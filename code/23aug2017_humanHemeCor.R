#+ echo=FALSE, message=FALSE, warning=FALSE
library(readr)
library(plotly)
library(heatmaply)
library(BuenColors)
library(ggplot2)
library(data.table)
library(GenomicRanges)

# Import human data
hemecounts <- data.matrix(fread("zcat < ../humandata/GSE74912_ATACseq_All_Counts.txt.gz"))
h <- hemecounts[,4:80]
 cc <- colnames(h)
 shortnames <- sapply(1:length(cc), function(s){
   if(s <= 23 | (s <= 74 & s >= 62 )){ strsplit(cc[s], split = "-")[[1]][2]
   } else { strsplit(cc[s], split = "-")[[1]][3]}
 })
noa <- gsub("A", "", shortnames)
nob <- gsub("B", "", noa)
noc <- gsub("C", "", nob)
#
n <- as.character(seq(1:15))
cell <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Mono", "missing", "CD4", "CD8", "NK", "missing", "Bcell", "CLP", "Ery")
names(cell) <- n
namess <- paste0(unname(cell[noc]))
colnames(h) <- namess

cellTypeDatMatHuman <- sapply(unique(namess), function(namee){
  rowSums(h[,c(which(namee== namess))])
})
remove(hemecounts)

humanPeaks <- fread("../humandata/peaks.bed")
humanPeaks$idx <- paste0("human", as.character(1:dim(humanPeaks)[1]))
human_g <-  makeGRangesFromDataFrame(humanPeaks, seqnames = "Chr", start.field = "Start",
                                     end.field = "End", keep.extra.columns = TRUE)


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

