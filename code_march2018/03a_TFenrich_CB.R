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

# Compile final results
resdf <- data.frame(genes = c(mcols(posg)$gene.name,mcols(negg)$gene.name),
                    data.matrix(motifMatches(matchm)),
                    stringsAsFactors = FALSE)


qdf <- read.table(vcFile, header = TRUE)
qdf$UnexplainedGene <- qdf$Unexplained > 90
qdf$PromoterGene <- qdf$Promoter > 90
qdf$DistalGene <- qdf$Distal > 90

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
distal$distal_stat <-  sign(log2(distal[,2]))* -1*log10(distal[,3])

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
promoter$promoter_stat <-  sign(log2(promoter[,1]))* -1*log10(promoter[,2])

all <- cbind(distal, promoter)
p1 <- ggplot(all, aes(x = distal_stat, promoter_stat, label = TF)) +
  geom_point(size = 0.5) + pretty_plot(fontsize = 6) + L_border() +
  labs(y = "Promoter Enrichment", x = "Distal Enrichment") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
cowplot::ggsave(p1, filename = "figures/scatter2C.pdf", width = 3, height = 3)
