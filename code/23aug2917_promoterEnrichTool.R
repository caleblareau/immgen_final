library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(BuenColors)

gtf <- "../data/mm10.refGenes.2016.1018.csv"
gdf <- fread(gtf)

featureMatch <- function(gdf, row, quiet = FALSE){
  
  # Define variables of interest
  feature <- as.character(row[1,1])
  far <- as.numeric(row[1,2])
  near <- as.numeric(row[1,3])
  motif <- as.character(row[1,4])
  
  # Quick check
  if(near > far){
    ttt <- near
    near <- far
    far <- ttt
  }
  
  # Process positive strand
  posdf <- gdf[gdf$strand == "+", c("chrom", "TSS", "gene.name")]
  posdf$start <- posdf$TSS - far
  posdf$end <- posdf$TSS - near
  posg <- makeGRangesFromDataFrame(posdf, keep.extra.columns = TRUE)
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, posg)
  pmatch <- !sapply(vmatchPattern(motif, seqs, fixed = FALSE)@ends, is.null)
  
  # Process negative strand 
  negdf <- gdf[gdf$strand == "-", c("chrom", "TSS", "gene.name")]
  negdf$start <- negdf$TSS + near
  negdf$end <- negdf$TSS + far
  negg <- makeGRangesFromDataFrame(negdf, keep.extra.columns = TRUE)
  seqs <- reverseComplement(getSeq(BSgenome.Mmusculus.UCSC.mm10, negg))
  nmatch <- !sapply(vmatchPattern(motif, seqs, fixed = FALSE)@ends, is.null)
  
  # Compile final results
  resdf <- data.frame(genes = c(mcols(posg)$gene.name,mcols(negg)$gene.name),
                      hit  = c(as.numeric(pmatch), as.numeric(nmatch)), stringsAsFactors = FALSE)
  names(resdf) <- c("genes", feature)
    
  # Some messages
  if(!quiet){
    message(feature)
    message(paste0("Positive strand match rate: ", as.character(round(sum(pmatch)/length(posg), 4))))
    message(paste0("Negative strand match rate: ", as.character(round(sum(nmatch)/length(negg), 4))))    
  }
  
  return(resdf[match(gdf$gene.name, resdf$genes),])
}

featuredf <- read.table("../data/promoterFeatures.txt", header = TRUE, stringsAsFactors = FALSE)
flo <- lapply(1:dim(featuredf)[1], function(i){
  featureMatch(gdf, featuredf[i,])
})

# Do logistic regressions
qdf <- read.table("../output/19aug_varianceComponents.txt", header = TRUE)

LRout <- sapply(1:length(flo), function(j){
  mdf <- merge(flo[[j]], qdf, by.x = "genes", by.y = "gene")
  
  # Fit logistic regression models
  m1 <- glm(mdf[,2] ~ mdf[,3], family = "binomial")
  m2 <- glm(mdf[,2] ~ mdf[,4], family = "binomial")
  m3 <- glm(mdf[,2] ~ mdf[,5], family = "binomial")
  
  # Get signed p-value
  sp1 <- (sign(summary(m1)$coefficient[2,1]) * -1 * log10(summary(m1)$coefficient[2,4]))
  sp2 <- (sign(summary(m2)$coefficient[2,1]) * -1 * log10(summary(m2)$coefficient[2,4]))
  sp3 <- (sign(summary(m3)$coefficient[2,1]) * -1 * log10(summary(m3)$coefficient[2,4]))
  c(sp1, sp2, sp3)
})

# Do Fisher test
FTout <- sapply(1:length(flo), function(j){
  mdf <- merge(flo[[j]], qdf, by.x = "genes", by.y = "gene")
  
  ft <- function(type){
     a <- mdf[,type] > 99.9 & mdf[,2] > 0
     b <-  mdf[,type] > 99.9 & mdf[,2] == 0
     c <- mdf[,type] < 99.9 & mdf[,2] > 0
     d <- mdf[,type] < 99.9 & mdf[,2] == 0
     ff <- fisher.test(matrix(c(sum(a),sum(b),sum(c),sum(d)), nrow = 2))
     unname(sign(log(ff$estimate))) * -1 * log10(ff$p.value)
  }
  
  # Get signed p-value
  sp1 <- ft(colnames(qdf)[c(1)])
  sp2 <- ft(colnames(qdf)[c(2)])
  sp3 <- ft(colnames(qdf)[c(3)])
  c(sp1, sp2, sp3)
})


LRout <- data.frame(LRout)
colnames(LRout) <- sapply(flo, function(ele){ colnames(ele)[2]})
LRout$VarianceComponent <- colnames(qdf)[c(1,2,3)]
plotdf <- reshape2::melt(LRout, id.vars = "VarianceComponent")
names(plotdf) <- c("VarianceComponent", "Motif", "Signed_Pvalue")

p1 <- ggplot() + pretty_plot() + 
  geom_bar(data = plotdf, aes(x=Motif, y=Signed_Pvalue, fill=VarianceComponent),
           stat = "identity", colour="black", position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("dodgerblue", "green3", "red")) + 
  theme( legend.position="bottom")+labs(colour = "Variance Component", fill = "Variance Component") +
  theme(axis.text.x = element_text(vjust = 0.25, angle = 90)) +labs(x = NULL, y = "Signed P-value")
print(p1)
ggsave(p1, filename = "../figures/ImmGen_VC_promoter.pdf")

