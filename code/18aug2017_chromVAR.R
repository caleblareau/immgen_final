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
peaks <- fread("../data/ImmGenATAC1219.peak.bed")
peaks <- makeGRangesFromDataFrame(peaks, seqnames = "V1", start.field = "V2", end.field = "V3")


counts <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))

data("mouse_pwms_v1") # also mouse_pwms_v1
motif_ix <- matchMotifs(mouse_pwms_v1, counts, genome = BSgenome.Mmusculus.UCSC.mm10)
counts <- addGCBias(counts,  genome = BSgenome.Mmusculus.UCSC.mm10)
dev <- computeDeviations(object = counts, annotations = motif_ix)

zscores <- assays(dev)[["z"]]
sd <- assays(dev)[["deviations"]]

write.table(zscores, file = "../chromVAR/zscoreData.csv",row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
write.table(data.frame(x = row.names(zscores)), file = "../chromVAR/zscoreRows.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = colnames(zscores)), file = "../chromVAR/zscoreCols.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)

write.table(sd, file = "../chromVAR/deviationsData.csv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = row.names(sd)), file = "../chromVAR/deviationsRows.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)
write.table(data.frame(x = colnames(sd)), file = "../chromVAR/deviationsCols.tsv",row.names = FALSE, col.names = FALSE, sep = ",",  quote = FALSE)

