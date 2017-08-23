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
remove(h)
humanPeaks <- fread("../humandata/peaks.bed")
humanPeaks$idx <- paste0("human", as.character(1:dim(humanPeaks)[1]))
human_g <-  makeGRangesFromDataFrame(humanPeaks, seqnames = "Chr", start.field = "Start",
                                     end.field = "End", keep.extra.columns = TRUE)

# Filter only human cell types
keepersH <- c("HSC", "MPP", "CLP", "Mono", "Bcell", "NK", "CD8", "CD4")

hCounts <- cellTypeDatMatHuman[,keepersH]
remove(cellTypeDatMatHuman)
colnames(hCounts) <- paste0("h", keepersH)

# Import ImmGen data
Icounts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
Ipeaks <- fread("../data/ImmGenATAC1219.peak.bed")
keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepidx <- Ipeaks[[4]] %in% keepers
Icounts <- Icounts[keepidx,]
Ipeaks <- Ipeaks[keepidx,]

# Filter for comparison Celltypes
mHSC <- c("LTHSC.34-.BM", "LTHSC.34+.BM")
mMPP <- c("MMP3.48+.BM", "MMP4.135+.BM")
mCLP <- c("proB.CLP.BM")
mMono <- c("Mo.6C-II-.Bl", "Mo.6C+II-.Bl")
mBcell <- c("B.Fo.Sp")
mNK <- c("NK.27-11b+.Sp", "NK.27+11b-.Sp", "NK.27+11b+.Sp")
mCD8 <- c("T.8.Nve.Sp")
mCD4 <- c("T.4.Nve.Sp")

allTypes <- list(mHSC, mMPP, mCLP, mMono, mBcell, mNK, mCD8, mCD4)

mCounts <- sapply(allTypes, function(typevec){
  if(length(typevec) > 1){
    rowSums(Icounts[,typevec])
  } else {
    Icounts[,typevec]
  }
})

colnames(mCounts) <- paste0("m", keepersH)
remove(Icounts)

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

immgen_LO_counts <- mCounts[unname(imap[overlapdf$immgenPeak]),]
human_LO_counts <- hCounts[unname(hmap[overlapdf$humanPeak]),]

mCPM <- sweep(immgen_LO_counts, 2, colSums(immgen_LO_counts), FUN="/") * 1000000
hCPM <- sweep(human_LO_counts, 2, colSums(human_LO_counts), FUN="/") * 1000000

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

keepMe <- rowVar(mCPM) > 3 & rowVar(hCPM) > 3
hCPM <- hCPM[keepMe,]
mCPM <- mCPM[keepMe,]

cormat <- cor(log2(mCPM+1), log2(hCPM+1))

heatplot <- ggplot(reshape2::melt(cormat), aes(x=Var2, y=Var1, fill = value)) +
  geom_tile() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black", family = "Helvetica", size = 10)) +
  coord_fixed(ratio=1) + scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  labs(fill='Correlation') + ggtitle("Filtered (soft) for variable peaks")

ggsave(heatplot, file = "../figures/mouse-human-heatmap2.pdf")


dim(hCPM)
