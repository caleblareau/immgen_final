library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(motifmatchr)
library(BuenColors)
library(chromVARmotifs)
library(magrittr)
library(SummarizedExperiment)


cpg <- read.table("../data/mm10_cpgIslandExt.txt")
expressedGenes <- read.table("../data/expressed_genes_SM.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
cpgg <- makeGRangesFromDataFrame(cpg, seqnames.field = "V2", start.field = "V3", end.field = "V4")
gtf <- "../data/mm10.refGenes.2016.1018.csv"
gdf <- fread(gtf)
gdf <- gdf[gdf[["gene.name"]]%in% expressedGenes, ]

motifPromoter <- function(gdf, far ){
  
  # Define promoter shape
  near <- 00
  
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
  inCpg <- 1:length(bothg) %in% queryHits(findOverlaps(bothg, cpgg))
  
  # Compile final results
  resdf <- data.frame(genes = c(mcols(posg)$gene.name,mcols(negg)$gene.name),
                      CPGi = inCpg,
                      stringsAsFactors = FALSE)
  return(resdf)
}
farr <- 1
resdf <- motifPromoter(gdf, far = farr)

qdf <- read.table("../output/RNA_varianceComponents.txt", header = TRUE)
qdf$UnexplainedGene <- qdf$Unexplained > 90
qdf$PromoterGene <- qdf$Promoter > 90
qdf$DistalGene <- qdf$Distal > 90

bdf <- merge(qdf, resdf, by.x = "gene", by.y = "genes")

df <- data.frame(
  hasCPGi = c(sum(bdf[,"DistalGene"] & bdf[,"CPGi"]), sum(bdf[,"PromoterGene"] & bdf[,"CPGi"]), sum(bdf[,"UnexplainedGene"] & bdf[,"CPGi"]),
              sum((!bdf[,"DistalGene"] & !bdf[,"PromoterGene"] & !bdf[,"UnexplainedGene"])& bdf[,"CPGi"])),
  noCPGi = c(sum(bdf[,"DistalGene"] & !bdf[,"CPGi"]), sum(bdf[,"PromoterGene"] & !bdf[,"CPGi"]), sum(bdf[,"UnexplainedGene"] & !bdf[,"CPGi"]),
             sum((!bdf[,"DistalGene"] & !bdf[,"PromoterGene"] & !bdf[,"UnexplainedGene"])& !bdf[,"CPGi"]))
)
rownames(df) <- c("Distal", "Promoter", "Unexplained", "Mixed")
df
plotdf <- data.frame(CpGiFreq = c(df[,1] / (df[,1] + df[,2])))
plotdf$GeneType<- c("Distal", "Promoter", "Unexplained", "Mixed")

p1 <- ggplot() + pretty_plot() + 
  geom_bar(data = plotdf, aes(x=GeneType, y=CpGiFreq*100, fill=GeneType),stat = "identity", color = "black") +
  scale_fill_manual(values=c("dodgerblue", "black", "green3", "red")) + 
  theme( legend.position="bottom")+labs(colour = "Variance Component", fill = "Variance Component") +
  theme(axis.text.x = element_text(vjust = 0.25, angle = 90)) +
  labs(x = "Gene Type", y = "% TSS Overlapping CpG Island")

p1
ggsave(p1,  file = paste0("../output/cpgi_overlap.pdf"))

