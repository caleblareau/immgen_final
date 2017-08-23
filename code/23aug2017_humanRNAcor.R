#+ echo=FALSE, message=FALSE, warning=FALSE
library(BuenColors)
library(ggplot2)
library(data.table)
library(GenomicRanges)

keepersH <- c("HSC", "MPP", "CLP", "Mono", "Bcell", "NK", "CD8", "CD4")

# Import human data
rnaraw <- fread("zcat < ../humandata/GSE74246_RNAseq_All_Counts.txt.gz")
genes <- as.character(data.frame(rnaraw[,1])[,1])
rnadat <- rnaraw[,-1]
types <- unique(unlist(strsplit(colnames(rnadat), split = "-"))[c(FALSE, TRUE)])[1:13]
alltypes <- unlist(strsplit(colnames(rnadat), split = "-"))[c(FALSE, TRUE)]
rnadat <- data.matrix(data.frame(rnaraw[,-1]))

collapsedMat <- sapply(types, function(type){
  rowSums(rnadat[,which(type == alltypes)])
})


log2norm <- log2(sweep(rnadat, 2, (colSums(rnadat) / 1000000), "/") + 1)


# Import ImmGen data
Icounts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
Ipeaks <- fread("../data/ImmGenATAC1219.peak.bed")
keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepidx <- Ipeaks[[4]] %in% keepers
Icounts <- Icounts[keepidx,]
Ipeaks <- Ipeaks[keepidx,]

# Add liftover
lo <- fread("../liftover/ImmGen_hg19peaks.bed")
lo_g <- makeGRangesFromDataFrame(lo[lo[["V4"]] %in% keepers,], seqnames = "V1", start.field = "V2",
                                 end.field = "V3", keep.extra.columns = TRUE)
  
# Intersect human and liftover
ov <- findOverlaps(lo_g, human_g)
overlapdf <- data.frame(immgenPeak = as.character(mcols(lo_g)[queryHits(ov),1]),
                        humanPeak = as.character(mcols(human_g)[subjectHits(ov),1]), stringsAsFactors = FALSE)

overlapdf <- overlapdf[!duplicated(overlapdf$immgenPeak),]
overlapdf <- overlapdf[!duplicated(overlapdf$humanPeak),]
#write.table(overlapdf, file = "../liftover/overlapPeaks.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

# No do the mapping
imap <- 1:dim(Ipeaks)[1]
names(imap) <- as.character(Ipeaks[[4]])
hmap <- 1:length(human_g)
names(hmap) <- as.character(humanPeaks$idx )

immgen_LO_counts <- Icounts[unname(imap[overlapdf$immgenPeak]),]
human_LO_counts <- cellTypeDatMatHuman[unname(hmap[overlapdf$humanPeak]),]

cormat <- cor(log2(immgen_LO_counts+1), log2(human_LO_counts+1))
heatmaply(cormat, colors = jdb_palette("brewer_spectra"))

