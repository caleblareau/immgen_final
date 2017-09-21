library(data.table)
library(GenomicRanges)
library(diffloop)
library(magrittr)
library(Matrix)

keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
lo <- fread("../liftover/ImmGen_hg19peaks.bed")
lo <- lo[lo[["V4"]] %in% keepers,]
lo$start <- lo$V2 - 125
lo$end <- lo$V2 + 125

lo_g <- makeGRangesFromDataFrame(lo, seqnames = "V1", start.field = "start",
                                 end.field = "end", keep.extra.columns = TRUE)

makeHitMat <- function(file){
  dt <- data.frame(fread(paste0("zcat < ", file), header = FALSE))
  dt <- dt[complete.cases(dt),]
  gwas_g <- makeGRangesFromDataFrame(dt, seqnames = "V1", start.field = "V2", end.field = "V3")
  as.numeric(1:length(lo_g) %in% queryHits(findOverlaps(lo_g, gwas_g)))
}

files <- list.files("../gwasRaw", full.names = TRUE)
names <- gsub("../gwasRaw/", "", files) %>% gsub(pattern = "_pruned_rsq_0.8_expanded_rsq_0.8.bed.gz", replacement = "")

hitMat250 <- sapply(files, makeHitMat)
checkdf <- data.frame(colnames(hitMat250), names)
tail(checkdf)

colnames(hitMat250) <- names
hitMat250 <- Matrix::Matrix(hitMat250[,colSums(hitMat250) >= 50], sparse = TRUE)
saveRDS(hitMat250, "../gwas_filtered/hitMat250_filt50.rds")


write.table(lo[["V4"]], file = "../gwas_filtered/ImmGenPeaksKeep.txt", sep = "\t",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

