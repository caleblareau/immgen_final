library(data.table)
library(motifmatchr)
library(chromVAR)
library(chromVARmotifs)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
"%ni%" <- Negate("%in%")

# Import ATAC and filter peaks
keepers <- data.table::fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepers <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))
bed <- data.frame(data.table::fread("zcat < ../data/ImmGenATAC1219.peak.bed.gz", header = FALSE, sep = "\t"))
gr <- makeGRangesFromDataFrame(bed, seqnames.field = "V1", start.field = "V2", end.field = "V3")
gr <- gr[keepers]
saveRDS(gr, "immgenPeaks.gr.rds")

# Match motifs
m <- matchMotifs(mouse_pwms_v2, gr, BSgenome.Mmusculus.UCSC.mm10, out = "scores", p.cutoff = 0.99)
