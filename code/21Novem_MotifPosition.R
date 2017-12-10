library(motifmatchr)
library(chromVARmotifs)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
library(magrittr)

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

peaksdf <- fread(paste0("zcat <", "../data/ImmGenATAC1219.peak.bed.gz"))
keepers <- fread(paste0("zcat < ", "../data/immgen_good_peaks.txt.gz"), header = FALSE)[["V1"]]
peaksdf <- data.frame(peaksdf[peaksdf$V4 %in% keepers,])
dim(peaksdf)

# Match motifs
data("mouse_pwms_v2")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
motif_ix <- motif_ix <- matchMotifs(mouse_pwms_v2, peaks, genome = "BSgenome.Mmusculus.UCSC.mm10", out = "position")

lapply(1:length(motif_ix), function(i){
  df <- data.frame(motif_ix[[i]])
  df$motif <- names(motif_ix)[i]
  df
}) %>% rbindlist() %>% as.data.frame() -> allMatchesDF

write.table(allMatchesDF, file = "../output/allMatchesPositions_HY21November.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

