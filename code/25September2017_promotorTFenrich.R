library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(motifmatchr)
library(BuenColors)
library(chromVARmotifs)
library(magrittr)
data("mouse_pwms_v1")

gtf <- "../data/mm10.refGenes.2016.1018.csv"
gdf <- fread(gtf)

motifPromoter <- function(gdf, far = 1000){
  
  # Define promoter shape
  near <- 0

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
  matchm <- motifmatchr::matchMotifs(mouse_pwms_v1, bothg, genome = BSgenome.Mmusculus.UCSC.mm10)
  
  # Compile final results
  resdf <- data.frame(genes = c(mcols(posg)$gene.name,mcols(negg)$gene.name),
                      data.matrix(motifMatches(matchm)),
                      stringsAsFactors = FALSE)
  return(resdf)
}
resdf <- motifPromoter(gdf)

qdf <- read.table("../output/19aug_varianceComponents.txt", header = TRUE)
qdf$UnexplainedGene <- qdf$Unexplained > 90
qdf$PromoterGene <- qdf$Promoter > 90
qdf$DistalGene <- qdf$Distal > 90

bdf <- merge(qdf, resdf, by.x = "gene", by.y = "genes")

sout <- sapply(names(mouse_pwms_v1), function(tf){
  OR_P <- function(type){
    tt <- bdf[,c(tf, type)]
    ft <- fisher.test(matrix(c( sum(tt[,1] & tt[,2]), sum(!tt[,1] & tt[,2]), sum(tt[,1] & !tt[,2]), sum(!tt[,1] & !tt[,2])),ncol =2))
    return(c(ft$estimate, ft$p.value))
  }
  c(OR_P("UnexplainedGene"), OR_P("PromoterGene"), OR_P("DistalGene"))
}) %>% t()
colnames(sout) <- c("UnexplainedOR", "UnexplainedP", "PromoterOR", "PromoterP", "DistalOR", "DistalP")

head(sout[order(sout[,"DistalP"]),], 20)
head(sout[order(sout[,"PromoterP"]),], 20)

p1 <- ggplot() + pretty_plot() + 
  geom_bar(data = plotdf, aes(x=Motif, y=Signed_Pvalue, fill=VarianceComponent),
           stat = "identity", colour="black", position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("dodgerblue", "green3", "red")) + 
  theme( legend.position="bottom")+labs(colour = "Variance Component", fill = "Variance Component") +
  theme(axis.text.x = element_text(vjust = 0.25, angle = 90)) +labs(x = NULL, y = "Signed P-value")
print(p1)
ggsave(p1, filename = "../figures/ImmGen_VC_promoter.pdf")

