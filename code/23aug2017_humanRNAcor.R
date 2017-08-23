#+ echo=FALSE, message=FALSE, warning=FALSE
library(BuenColors)
library(ggplot2)
library(data.table)
library(GenomicRanges)

keepersH <- c("HSC", "MPP", "CLP", "Mono", "Bcell", "NK", "CD8", "CD4")

# Import human data
rnaraw <- fread("zcat < ../humandata/GSE74246_RNAseq_All_Counts.txt.gz")
hGenes <- as.character(data.frame(rnaraw[,1])[,1])
rnadat <- rnaraw[,-1]
types <- unique(unlist(strsplit(colnames(rnadat), split = "-"))[c(FALSE, TRUE)])[1:13]
alltypes <- unlist(strsplit(colnames(rnadat), split = "-"))[c(FALSE, TRUE)]
rnadat <- data.matrix(data.frame(rnaraw[,-1]))

hCounts <- sapply(c("HSC", "MPP", "CLP", "Mono", "Bcell", "NKcell", "CD8Tcell", "CD4Tcell"), function(type){
  rowSums(rnadat[,which(type == alltypes)])
})
remove(rnadat)
colnames(hCounts) <- paste0("h", keepersH)


# Mouse defintions
mHSC <- c("LTHSC.34-.BM", "LTHSC.34+.BM")
mMPP <- c("MMP3.48+.BM", "MMP4.135+.BM")
mCLP <- c("proB.CLP.BM")
mMono <- c("Mo.6C-II-.Bl", "Mo.6C+II-.Bl")
mBcell <- c("B.Fo.Sp")
mNK <- c("NK.27-11b+.Sp", "NK.27+11b-.Sp", "NK.27+11b+.Sp")
mCD8 <- c("T.8.Nve.Sp")
mCD4 <- c("T.4.Nve.Sp")

allTypes <- list(mHSC, mMPP, mCLP, mMono, mBcell, mNK, mCD8, mCD4)

# Import ImmGen data
Icounts <- data.matrix(fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
mGenes <- read.table("../data/rnaRowLabels.txt")[,1]

mCounts <- sapply(allTypes, function(typevec){
  if(length(typevec) > 1){
    rowSums(Icounts[,typevec])
  } else {
    Icounts[,typevec]
  }
})
remove(Icounts)
colnames(mCounts) <- paste0("m", keepersH)

# Add liftover gene IDs
lo_dt <- fread("../liftover/HOM_MouseHumanSequence.rpt")
transdf <- data.frame(mouse = lo_dt[["Symbol"]][c(TRUE,FALSE)], human = lo_dt[["Symbol"]][c(FALSE,TRUE)], stringsAsFactors = FALSE)
mouseToHuman <- transdf$human
names(mouseToHuman) <- transdf$mouse
humanToMouse <- transdf$mouse
names(humanToMouse) <- transdf$human

# Make a consensus set of genes
removeNA <- function(x) x[!is.na(x)]
hmhGenes  <- mouseToHuman[removeNA(unname(humanToMouse[hGenes]))]
hmGenes <- removeNA(unname(mouseToHuman[mGenes]))

keepHuman <- sort(intersect(hmGenes, hmhGenes))
keepMouse <- humanToMouse[keepHuman]

masterTranslateDF <- data.frame(
  humanGene = keepHuman, 
  mouseGene = keepMouse, 
  humanIdx = match(keepHuman, hGenes),
  mouseIdx = match(keepMouse, mGenes),
  stringsAsFactors = FALSE
)
masterTranslateDF <- masterTranslateDF[complete.cases(masterTranslateDF), ]
dim(masterTranslateDF)

immgen_LO_counts <- mCounts[masterTranslateDF$mouseIdx,]
human_LO_counts <- hCounts[masterTranslateDF$humanIdx,]

mCPM <- sweep(immgen_LO_counts, 2, colSums(immgen_LO_counts), FUN="/") * 1000000
hCPM <- sweep(human_LO_counts, 2, colSums(human_LO_counts), FUN="/") * 1000000

cormat <- cor(log2(mCPM+1), log2(hCPM+1), use = "com")

heatplot <- ggplot(reshape2::melt(cormat), aes(x=Var2, y=Var1, fill = value)) +
  geom_tile() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black", family = "Helvetica", size = 10)) +
  coord_fixed(ratio=1) + scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  labs(fill='Correlation') + ggtitle("No filtering for variable peaks")
print(heatplot)


ggsave(heatplot, file = "../figures/mouse-human-RNA-heatmap0.pdf")



