library(irlba)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(cowplot)
library(data.table)
library(GenomicRanges)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)

"%ni%" <- Negate("%in%")

blackout <- c("Ep.MEChi.Th", "FRC.CD140a+.Madcam-.CD35-.SLN", "LEC.SLN", "BEC.SLN", "IAP.SLN")

keepers <- fread(paste0("zcat < ", "../data/immgen_good_peaks.txt.gz"), header = FALSE)[[1]]
keepers <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))

atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
atac <- atac[keepers,colnames(atac) %ni% blackout]

peaksdf <- fread(paste0("zcat < ", "../data/ImmGenATAC1219.peak.bed.gz"))
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <- SummarizedExperiment(assays = list(counts = atac),
                               rowData = peaks[keepers,], 
                               colData = DataFrame(names = colnames(atac)))
counts <- addGCBias(counts,  genome = BSgenome.Mmusculus.UCSC.mm10)

bg <- getBackgroundPeaks(counts, niterations = 200)
write.table(bg, file = "../output/bgpeaks_200.csv", row.names = FALSE, col.names = FALSE, sep=",", quote = FALSE)
