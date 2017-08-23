library(data.table)

ovPeaks <- fread("../liftover/overlapPeaks.txt", header = FALSE)
keepPeaks <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
loPeaks <- as.numeric(gsub("ImmGenATAC1219.peak_", "", ovPeaks[[1]]))
peakPvals <- data.matrix(fread("zcat < ../data/0726.ATAC.population.log10Pval.BH.csv.gz")[,V1:=NULL])
loPeaksPvals <- peakPvals[loPeaks,]
peakPvals <- peakPvals[as.numeric(gsub("ImmGenATAC1219.peak_", "", keepPeaks)),]
all(ovPeaks[[1]] %in% keepPeaks)


getLiftoverRatio <- function(threshold){
  num <- colSums(loPeaksPvals > threshold)
  denom <- colSums(peakPvals > threshold)
  return(num/denom)
}


odf <- data.frame(
  sample = colnames(peakPvals),
  FDR05 = getLiftoverRatio(-1*log10(0.05)),
  FDR01 = getLiftoverRatio(-1*log10(0.01)),
  FDR001 = getLiftoverRatio(-1*log10(0.001))
)

write.table(odf, file = "../liftover/populationPercentages/liftoverEfficacy_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE, 
            col.names = TRUE)


library(data.table)

ovPeaks <- fread("../liftover/overlapPeaks.txt", header = FALSE)
keepPeaks <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepPeaks <- as.numeric(gsub("ImmGenATAC1219.peak_", "", keepPeaks))
loPeaks <- as.numeric(gsub("ImmGenATAC1219.peak_", "", ovPeaks[[1]]))
peakPvals <- data.matrix(fread("zcat < ../data/0726.ATAC.population.log10Pval.BH.csv.gz")[,V1:=NULL])

peakAnno <- fread("../data/peakAnnoVector.txt")[[1]]

makeOutTable <- function(annotation){
  
  inAnnotation <- peakAnno == annotation
  
  xloPeaksPvals <- peakPvals[1:dim(peakPvals)[1] %in% loPeaks & inAnnotation,]
  xpeakPvals <- peakPvals[1:dim(peakPvals)[1] %in% keepPeaks & inAnnotation,]

  getLiftoverRatio <- function(threshold){
    num <- colSums(xloPeaksPvals > threshold)
    denom <- colSums(xpeakPvals > threshold)
    return(num/denom)
  }
  
  odf <- data.frame(
  sample = colnames(peakPvals),
  FDR05 = getLiftoverRatio(-1*log10(0.05)),
  FDR01 = getLiftoverRatio(-1*log10(0.01)),
  FDR001 = getLiftoverRatio(-1*log10(0.001))
)

write.table(odf, file = paste0("../liftover/populationPercentages/liftoverEfficacy_",annotation,".tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}
makeOutTable("TSS")
makeOutTable("Genebody")
makeOutTable("outside")
