library(data.table)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(chromVARmotifs)
library(GenomicRanges)
library(BuenColors)
library(diffloop)
library(magrittr)

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

counts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
peaksdf <- fread("../data/ImmGenATAC1219.peak.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
keepers <- fread("../gwas_filtered/ImmGenPeaksKeep.txt", header = FALSE)[["V1"]]

counts <- SummarizedExperiment(assays = list(counts = counts[peaksdf[["V4"]] %in% keepers,]),
                               rowData = peaks[peaksdf[["V4"]] %in% keepers,], 
                               colData = DataFrame(names = colnames(counts)))

counts <- addGCBias(counts,  genome = BSgenome.Mmusculus.UCSC.mm10)

# Build GWAS hits matrix
keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
lo <- fread("../liftover/ImmGen_hg19peaks.bed")
lo <- lo[lo[["V4"]] %in% keepers,]
lo$start <- lo$V2 - 125
lo$end <- lo$V2 + 125

lo_g <- makeGRangesFromDataFrame(lo, seqnames = "V1", start.field = "start",
                                 end.field = "end", keep.extra.columns = TRUE)


bimdf <- fread(paste0("zcat <", "../data/snpcoords.txt.gz"))
makeHitMat <- function(file){
  df <- read.table(file, header = FALSE)
  dt <- merge(df, bimdf, by = "V1")
  gwas_g <- makeGRangesFromDataFrame(dt, seqnames = "V2", start.field = "V3", end.field = "V3")
  as.numeric(1:length(lo_g) %in% queryHits(findOverlaps(rmchr(lo_g), gwas_g)))
}

files <- list.files("../gwasRaw", full.names = TRUE)
gwas_ix <- sapply(files, makeHitMat)
colnames(gwas_ix) <- gsub(".tags", "", gsub("../gwasRaw/", "", files))

colSums(gwas_ix)
dev <- computeDeviations(object = counts, annotations = gwas_ix)

zscores <- assays(dev)[["z"]]
zm <- reshape2::melt(zscores)
heatmaply::heatmaply(zscores, colors = jdb_palette("brewer_spectra"))
