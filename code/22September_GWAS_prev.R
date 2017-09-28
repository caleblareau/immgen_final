library(data.table)
library(GenomicRanges)
library(diffloop)
library(magrittr)
library(Matrix)

"%ni%" <- Negate("%in%")


keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
lo <- fread("../liftover/ImmGen_hg19peaks.bed")
lo <- lo[lo[["V4"]] %in% keepers,]
lo$start <- (lo$V2 + lo$V3)/2 - 125
lo$end <- (lo$V2 + lo$V3)/2 + 125

lo_g <- makeGRangesFromDataFrame(lo, seqnames = "V1", start.field = "start",
                                 end.field = "end", keep.extra.columns = TRUE)

humanPeaks <- fread("../humandata/peaks.bed")
humanPeaks$idx <- paste0("human", as.character(1:dim(humanPeaks)[1]))
human_g <-  makeGRangesFromDataFrame(humanPeaks, seqnames = "Chr", start.field = "Start",
                                     end.field = "End", keep.extra.columns = TRUE)

defined <- Matrix(data.matrix(fread("zcat < ../gwasRaw/hema_peaks_w_NONSIG_gwas_hits_above_thresh1e-05.txt.gz")))
ov <- findOverlaps(lo_g, human_g)

laout <- sapply(1:dim(lo)[1], function(i) {
  if( i %ni% queryHits(ov)){
    return(rep(0, dim(defined)[2]))
  } else if (length(which(i == queryHits(ov))) == 1){
    unname(defined[subjectHits(ov)[which(i == queryHits(ov))],])
  } else {
    unname(Matrix::colSums( defined[subjectHits(ov)[which(i == queryHits(ov))],]))
  }
}) %>% t() %>% Matrix()

colnames(laout) <- colnames(defined)
saveRDS(laout, "../gwas_filtered/22September_new.rds")

