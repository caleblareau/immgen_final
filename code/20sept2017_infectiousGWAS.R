library(data.table)
library(GenomicRanges)
library(diffloop)
library(magrittr)
library(Matrix)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BuenColors)

BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

# Import the peaks and translation component
keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
lo <- fread("../liftover/ImmGen_hg19peaks.bed")
lo <- lo[lo[["V4"]] %in% keepers,]
lo$start <- lo$V2 - 125
lo$end <- lo$V2 + 125

lo_g <- makeGRangesFromDataFrame(lo, seqnames = "V1", start.field = "start",
                                 end.field = "end", keep.extra.columns = TRUE)

# Deal with the infectious disease stuff
infraw <- fread("../data/infectious_sumSTATs.txt")
infraw <- infraw[!(infraw$scaffold == "chr6" & infraw$position > 29691116 & 33054976 < infraw$position ), ] # Remove HLA
makeHitMatInf <- function(trait){
  gwas_g <- makeGRangesFromDataFrame(data.frame(infraw[infraw$phenotype == trait,]),
                                     seqnames = "scaffold", start.field = "position", end.field = "position")
  as.numeric(1:length(lo_g) %in% queryHits(findOverlaps(lo_g, gwas_g)))
}
hitMatInf <- sapply( unique(infraw[["phenotype"]]), makeHitMatInf)

# Iport counts data
counts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
peaksdf <- fread("../data/ImmGenATAC1219.peak.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
keepers <- fread("../gwas_filtered/ImmGenPeaksKeep.txt", header = FALSE)[["V1"]]

counts <- SummarizedExperiment(assays = list(counts = counts[peaksdf[["V4"]] %in% keepers,]),
                               rowData = peaks[peaksdf[["V4"]] %in% keepers,], 
                               colData = DataFrame(names = colnames(counts)))

counts <- addGCBias(counts,  genome = BSgenome.Mmusculus.UCSC.mm10)
gwas_ix <- hitMatInf
devInf <- computeDeviations(object = counts, annotations = gwas_ix)

zscores <- assays(devInf)[["z"]]
zscores[zscores > 3] <- 3
zscores[zscores < -3] <- -3
heatmaply::heatmaply(zscores, colors = jdb_palette("brewer_spectra"))

