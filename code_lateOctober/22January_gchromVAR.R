library(gchromVAR)
library(data.table)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(chromVARmotifs)
library(GenomicRanges)
library(BuenColors)

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

# Build Counts Matrix
counts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
peaksdf <- fread("zcat < ../data/ImmGenATAC1219.peak.bed.gz")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
keepers <- fread("../gwas_filtered/ImmGenPeaksKeep.txt", header = FALSE)[["V1"]]

SE <- SummarizedExperiment(assays = list(counts = counts[peaksdf[["V4"]] %in% keepers,]),
                               rowData = peaks[peaksdf[["V4"]] %in% keepers,], 
                               colData = DataFrame(names = colnames(counts)))

SE <- addGCBias(SE,  genome = BSgenome.Mmusculus.UCSC.mm10)

# Build GWAS Matrix
gwasSE <- importBedScore(rowRanges(SE), list.files("../pics_bedscore/mm10/", full.names = TRUE))
assayNames(gwasSE) <- "matches"


# Run gchromVAR
wdev <- computeDeviations(SE, gwasSE)
widemat <- t(assays(wdev)[["z"]])
df <- melt(widemat)
head(df[order(df$value, decreasing = TRUE),], 20)

library(dplyr)
df %>% group_by(Var2) %>% top_n(n = 1, wt = value) %>% as.data.frame()

heatplot <- ggplot(reshape2::melt(widemat), aes(x=Var2, y=Var1)) +
                geom_tile(aes(fill = value),colour = "white") + theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.text = element_text(colour = "black", family = "Helvetica", size = 10)) +
                 coord_fixed(ratio=1) + scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                theme(axis.title.x=element_blank(), axis.title.y=element_blank())
cowplot::ggsave(heatplot, filename = "../figures/22January_gchromVAR.pdf", width = 30, height = 45)


