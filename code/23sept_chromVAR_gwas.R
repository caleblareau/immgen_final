#+ echo=FALSE, message=FALSE, warning=FALSE
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

counts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
peaksdf <- fread("../data/ImmGenATAC1219.peak.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
keepers <- fread("../gwas_filtered/ImmGenPeaksKeep.txt", header = FALSE)[["V1"]]

counts <- SummarizedExperiment(assays = list(counts = counts[peaksdf[["V4"]] %in% keepers,]),
                               rowData = peaks[peaksdf[["V4"]] %in% keepers,], 
                               colData = DataFrame(names = colnames(counts)))

counts <- addGCBias(counts,  genome = BSgenome.Mmusculus.UCSC.mm10)
gwas_ix <- readRDS("../gwas_filtered/22September_new.rds")
gwas_ix <- gwas_ix[, Matrix::colSums(gwas_ix) > 50]
inblood <- gsub("-", "_", read.table("../gwas_filtered/enrichedInBlood_corces.txt", header = FALSE)[,1])
dev <- computeDeviations(object = counts, annotations = gwas_ix)

zscores <- assays(dev)[["z"]]
#zscores[zscores > 3] <- 3
#zscores[zscores < -3] <- -3

neat <- c("Alzheimers_disease", "Celiac_disease", "Rheumatoid_arthritis", "Type_1_diabetes", "psoriasis_gwas_HGVRS1736_sorted",
          "Irritible_bowel_syndrome", "EUR.UC", "EUR.CD" )

zscores <- zscores[neat,]
heatmaply::heatmaply(zscores, colors = jdb_palette("brewer_spectra"))

odf <- data.frame(trait = rownames(assays(dev)[["z"]]),
                  blood_enriched = rownames(assays(dev)[["z"]]) %in% inblood,
                  assays(dev)[["z"]])

write.table(odf, file = "../gwas_filtered/20September_chromVAR_gwas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

zscores <- read.table("../gwas_filtered/20September_chromVAR_gwas.tsv", header = TRUE)
rownames(zscores) <- zscores[,1]
zscores <- zscores[,c(-1,-2)]
zscores[zscores > 3] <- 3
zscores[zscores < -3] <- -3

heatmaply::heatmaply(zscores, colors = jdb_palette("brewer_spectra"))

what <- "Multiple"
rownames(odf)[grep(what, rownames(odf))]
sort(odf[rownames(odf)[grep(what, rownames(odf))][1],])

hmm <- apply(odf, 1, function(r) names(sort(r, decreasing = TRUE)[3]))

