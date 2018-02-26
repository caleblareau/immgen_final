#+ echo=FALSE, message=FALSE, warning=FALSE
library(data.table)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(chromVARmotifs)
library(GenomicRanges)


# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

#' Analysis of ImmGen Consortium ATAC data using chromVAR
# Performed by Caleb Lareau, 18 August

#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
counts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
peaks <- fread("zcat < ../data/ImmGenATAC1219.peak.bed.gz")

keepers <- fread(paste0("zcat < ", "../data/immgen_good_peaks.txt.gz"), header = FALSE)[[1]]
keepers <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))

counts <- counts[keepers,]
peaks <- peaks[keepers,]

peaks <- makeGRangesFromDataFrame(peaks, seqnames = "V1", start.field = "V2", end.field = "V3")

counts <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
counts <- addGCBias(counts,  genome = BSgenome.Mmusculus.UCSC.mm10)
bg <- getBackgroundPeaks(counts)
#write.table(bg, file = "../chromVAR/18aug_backgroundPeaks.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
remove(bg)

# Motifs
data("mouse_pwms_v2") # also mouse_pwms_v1
motif_ix <- matchMotifs(mouse_pwms_v2, counts, genome = BSgenome.Mmusculus.UCSC.mm10)

# motifdf <- data.frame(data.matrix(assays(motif_ix)[["motifMatches"]]))
# write.table(motifdf, file = "../chromVAR/motifMatchData.csv",row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
# write.table(data.frame(x = colnames(motifdf)), file = "../chromVAR/motifMatchCols.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)


# Scores
dev <- computeDeviations(object = counts, annotations = motif_ix)
zscores <- assays(dev)[["z"]]
sd <- assays(dev)[["deviations"]]
variability <- computeVariability(dev)

write.table(variability, file = "../chromVAR/variabilityData.csv",row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
write.table(data.frame(x = row.names(variability)), file = "../chromVAR/variabilityRows.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = colnames(variability)), file = "../chromVAR/variabilityCols.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)

write.table(zscores, file = "../chromVAR/zscoreData.csv",row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
write.table(data.frame(x = row.names(zscores)), file = "../chromVAR/zscoreRows.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = colnames(zscores)), file = "../chromVAR/zscoreCols.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)

write.table(sd, file = "../chromVAR/deviationsData.csv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = row.names(sd)), file = "../chromVAR/deviationsRows.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = colnames(sd)), file = "../chromVAR/deviationsCols.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)

