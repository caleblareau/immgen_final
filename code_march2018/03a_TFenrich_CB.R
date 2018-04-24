library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(motifmatchr)
library(BuenColors)
library(chromVARmotifs)
library(magrittr)
library(SummarizedExperiment)
library(broom)
library(dplyr)

"%ni%" <- Negate("%in%")

data("mouse_pwms_v2")
vcFile <- "output/varianceComponentsOutput.txt"
gtf <- "../data/mm10.refGenes.2016.1018.csv"
gdf <- fread(gtf)


# Define promoter shape
near <- 0
far <- 1000

# Process positive strand
posdf <- gdf[gdf$strand == "+", c("chrom", "TSS", "gene.name")]
posdf$start <- posdf$TSS - far
posdf$end <- posdf$TSS - near
posg <- makeGRangesFromDataFrame(posdf, keep.extra.columns = TRUE)

# Process negative strand 
negdf <- gdf[gdf$strand == "-", c("chrom", "TSS", "gene.name")]
negdf$start <- negdf$TSS + near
negdf$end <- negdf$TSS + far
negg <- makeGRangesFromDataFrame(negdf, keep.extra.columns = TRUE)

bothg <- c(posg, negg)
matchm <- motifmatchr::matchMotifs(mouse_pwms_v2, bothg, genome = BSgenome.Mmusculus.UCSC.mm10)
SE <- SummarizedExperiment(
  rowData = bothg
  )
gc <- chromVAR::addGCBias( SE, genome = BSgenome.Mmusculus.UCSC.mm10)

# Compile final results
resdf <- data.frame(genes = c(mcols(posg)$gene.name,mcols(negg)$gene.name),
                    data.matrix(motifMatches(matchm)),
                    stringsAsFactors = FALSE)
g = mcols(gc)$bias

qdf <- read.table(vcFile, header = TRUE)
qdf$UnexplainedGene <- qdf$Unexplained > 95
qdf$PromoterGene <- qdf$Promoter > 95
qdf$DistalGene <- qdf$Distal > 95
genes <- as.character(mcols(gc)$gene.name)

dg <- as.character(vco$gene)[qdf$DistalGene]
pg <- as.character(vco$gene)[qdf$PromoterGene] 
ug <- as.character(vco$gene)[qdf$UnexplainedGene ] 

gcdf <- data.frame(
  GC = c(g[genes %in% dg], g[genes %in% pg], g[genes %in% ug]),
  Type = c(rep("Distal", sum(genes %in% dg)), rep("Promoter", sum(genes %in% pg)), rep("Unexplained",sum(genes%in% ug)))
)

p3 <- ggplot(gcdf, aes(x = Type, y = GC, color = Type)) + geom_boxplot() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "", y = "Promoter GC Content",  color = "") +
  scale_color_manual(values = c("dodgerblue3", "green3", "red"))

cowplot::ggsave(cowplot::plot_grid(p3, nrow = 1), filename = "figures/promoterGCcontent.pdf", width = 4, height = 2.5)


genes <- qdf[qdf$DistalGene,"gene"]
lapply(2:dim(resdf)[2], function(i){
  a <- sum(resdf[resdf$genes %in% genes,i])
  b <- sum(!resdf[resdf$genes %in% genes,i])
  c <-  sum(resdf[resdf$genes %ni% genes,i])
  d <- sum(!resdf[resdf$genes %ni% genes,i])
  ft <- tidy(fisher.test(matrix(c(a,b,c,d),
                                nrow =2)))[,c(1,2)]
  data.frame(colnames(resdf)[i], ft)
}) %>% rbindlist() %>% data.frame() -> distal
colnames(distal) <- c("TF", "distal_estimate", "distal_pvalue")
distal$distal_stat <-  round(sign(log2(distal[,2]))* -1*log10(distal[,3]), 1)

genes <- qdf[qdf$PromoterGene,"gene"]
lapply(2:dim(resdf)[2], function(i){
  a <- sum(resdf[resdf$genes %in% genes,i])
  b <- sum(!resdf[resdf$genes %in% genes,i])
  c <-  sum(resdf[resdf$genes %ni% genes,i])
  d <- sum(!resdf[resdf$genes %ni% genes,i])
  ft <- tidy(fisher.test(matrix(c(a,b,c,d),
                                nrow =2)))[,c(1,2)]
  data.frame(ft)
}) %>% rbindlist() %>% data.frame() -> promoter
colnames(promoter) <- c("promoter_estimate", "promoter_pvalue")
promoter$promoter_stat <-  round(sign(log2(promoter[,1]))* -1*log10(promoter[,2]), 1)

all <- cbind(distal, promoter)
p1 <- ggplot(all, aes(x = distal_stat, promoter_stat, label = TF)) +
  geom_point(size = 0.5) + pretty_plot(fontsize = 6) + L_border() +
  labs(y = "Promoter Enrichment", x = "Distal Enrichment") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
cowplot::ggsave(p1, filename = "figures/scatter2C.pdf", width = 3, height = 3)

gcdf 

write.table(all[,c(1,4,7)], file = paste0("output/varianceComponentMotifs_april20.tsv"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
