library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(motifmatchr)
library(BuenColors)
library(chromVARmotifs)
library(magrittr)
library(SummarizedExperiment)

data("mouse_pwms_v1")

gtf <- "../data/mm10.refGenes.2016.1018.csv"
gdf <- fread(gtf)

doAll <- function(farr){
  motifPromoter <- function(gdf, far ){
    
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
  resdf <- motifPromoter(gdf, far = farr)
  
  qdf <- read.table("../output/19aug_varianceComponents.txt", header = TRUE)
  qdf$UnexplainedGene <- qdf$Unexplained > 90
  qdf$PromoterGene <- qdf$Promoter > 90
  qdf$DistalGene <- qdf$Distal > 90
  
  bdf <- merge(qdf, resdf, by.x = "gene", by.y = "genes")
  bdf2 <- bdf[bdf$PromoterGene | bdf$DistalGene,]
  
  enrichOut <- sapply(names(mouse_pwms_v1), function(tf){
    OR_P <- function(type){
      tt <- bdf[,c(tf, type)]
      ft <- fisher.test(matrix(c( sum(tt[,1] & tt[,2]), sum(!tt[,1] & tt[,2]), sum(tt[,1] & !tt[,2]), sum(!tt[,1] & !tt[,2])),ncol =2))
      return(c(ft$estimate, ft$p.value))
    }
    c(OR_P("DistalGene"))
  }) %>% t()
  
  plotdf <- data.frame(stat = sign(log2(enrichOut[,1]))* -1*log10(enrichOut[,2]),
                       tf = rownames(enrichOut), sign = sign(log2(enrichOut[,1])))
  plotdf$sign <- ifelse(plotdf$sign == 1, "Distal", "Promoter")
  plotdf <- plotdf[order(plotdf$stat, decreasing = TRUE),]
  plotdf$rank <- 1:dim(plotdf)[1]
  plotdf <- plotdf[plotdf$sign != 0,]
  
  enrichOutPermuted <- sapply(names(mouse_pwms_v1), function(tf){
    OR_P <- function(type){
      tt <- data.frame(bdf[,c(tf)], sample(bdf[,c(tf)]))
      ft <- fisher.test(matrix(c( sum(tt[,1] & tt[,2]), sum(!tt[,1] & tt[,2]), sum(tt[,1] & !tt[,2]), sum(!tt[,1] & !tt[,2])),ncol =2))
      return(c(ft$estimate, ft$p.value))
    }
    c(OR_P("DistalGene"))
  }) %>% t()
  
  permutedf <- data.frame(stat = sign(log2(enrichOutPermuted[,1]))* -1*log10(enrichOutPermuted[,2]),
                          tf = rownames(enrichOutPermuted))
  permutedf$sign <- "Permuted"
  permutedf <- permutedf[order(permutedf$stat, decreasing = TRUE),]
  permutedf$rank <- 1:dim(permutedf)[1]
  permutedf <- permutedf[permutedf$sign != 0,]
  
  write.table(plotdf, file = paste0("../output/tfenrich/varianceComponent_promoterTFenrich.", as.character(farr), "bp.results.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  totaldf <- rbind(plotdf,permutedf)
  totaldf$sign <- factor(totaldf$sign, levels = c("Distal", "Promoter", "Permuted"))
  p1 <- ggplot() + pretty_plot() + 
    geom_bar(data =totaldf, aes(x=rank, y=stat, fill=sign),stat = "identity") +
    scale_fill_manual(values=c("dodgerblue", "green3", "grey")) + 
    theme( legend.position="bottom")+labs(colour = "Variance Component", fill = "Variance Component") +
    theme(axis.text.x = element_text(vjust = 0.25, angle = 90)) +
    labs(x = "Rank Sorted Transcription Factor Motif", y = "Signed -log10p")
  ggsave(p1,  file = paste0("../output/tfenrich/varianceComponent_promoterTFenrich.", as.character(farr), "bp.pdf"))
}


doKmer <- function(farr){
  motifPromoter <- function(gdf, far){
    
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
    matchm <- chromVAR::matchKmers(6, bothg, genome = BSgenome.Mmusculus.UCSC.mm10)
    
    # Compile final results
    resdf <- data.frame(genes = c(mcols(posg)$gene.name,mcols(negg)$gene.name),
                        data.matrix(assays(matchm)[["matches"]]),
                        stringsAsFactors = FALSE)
    return(resdf)
  }
  resdf <- motifPromoter(gdf, far = farr)
  
  qdf <- read.table("../output/19aug_varianceComponents.txt", header = TRUE)
  qdf$UnexplainedGene <- qdf$Unexplained > 90
  qdf$PromoterGene <- qdf$Promoter > 90
  qdf$DistalGene <- qdf$Distal > 90
  
  bdf <- merge(qdf, resdf, by.x = "gene", by.y = "genes")
  bdf2 <- bdf[bdf$PromoterGene | bdf$DistalGene,]
  
  kmers <- colnames( bdf2[,8:dim(bdf2)[2]])
  enrichOut <- sapply(kmers, function(tf){
    OR_P <- function(type){
      tt <- bdf[,c(tf, type)]
      ft <- fisher.test(matrix(c( sum(tt[,1] & tt[,2]), sum(!tt[,1] & tt[,2]), sum(tt[,1] & !tt[,2]), sum(!tt[,1] & !tt[,2])),ncol =2))
      return(c(ft$estimate, ft$p.value))
    }
    c(OR_P("DistalGene"))
  }) %>% t()
  
  plotdf <- data.frame(stat = sign(log2(enrichOut[,1]))* -1*log10(enrichOut[,2]),
                       tf = kmers, sign = sign(log2(enrichOut[,1])))
  plotdf$sign <- ifelse(plotdf$sign == 1, "Distal", "Promoter")
  plotdf <- plotdf[order(plotdf$stat, decreasing = TRUE),]
  plotdf$rank <- 1:dim(plotdf)[1]
  plotdf <- plotdf[plotdf$sign != 0,]
  
  enrichOutPermuted <- sapply(kmers, function(tf){
    OR_P <- function(type){
      tt <- data.frame(bdf[,c(tf)], sample(bdf[,c(tf)]))
      ft <- fisher.test(matrix(c( sum(tt[,1] & tt[,2]), sum(!tt[,1] & tt[,2]), sum(tt[,1] & !tt[,2]), sum(!tt[,1] & !tt[,2])),ncol =2))
      return(c(ft$estimate, ft$p.value))
    }
    c(OR_P("DistalGene"))
  }) %>% t()
  
  permutedf <- data.frame(stat = sign(log2(enrichOutPermuted[,1]))* -1*log10(enrichOutPermuted[,2]),
                          tf = rownames(enrichOutPermuted))
  permutedf$sign <- "Permuted"
  permutedf <- permutedf[order(permutedf$stat, decreasing = TRUE),]
  permutedf$rank <- 1:dim(permutedf)[1]
  permutedf <- permutedf[permutedf$sign != 0,]
  
  write.table(plotdf, file = paste0("../output/tfenrich/varianceComponent_promoterKMERenrich.", as.character(farr), "bp.results.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  totaldf <- rbind(plotdf,permutedf)
  totaldf$sign <- factor(totaldf$sign, levels = c("Distal", "Promoter", "Permuted"))
  p1 <- ggplot() + pretty_plot() + 
    geom_bar(data =totaldf, aes(x=rank, y=stat, fill=sign),stat = "identity") +
    scale_fill_manual(values=c("dodgerblue", "green3", "grey")) + 
    theme( legend.position="bottom")+labs(colour = "Variance Component", fill = "Variance Component") +
    theme(axis.text.x = element_text(vjust = 0.25, angle = 90)) +
    labs(x = "Rank Sorted K-mer", y = "Signed -log10p")
  ggsave(p1,  file = paste0("../output/tfenrich/varianceComponent_promoterKMERenrich.", as.character(farr), "bp.pdf"))
}

doAll(farr = 50)
doAll(farr = 100)
doAll(farr = 1000)

doKmer(farr = 50)
doKmer(farr = 100)
doKmer(farr = 1000)