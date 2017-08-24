library(irlba)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(cowplot)
library(data.table)
"%ni%" <- Negate("%in%")

blackout <- c("Ep.MEChi.Th", "FRC.CD140a+.Madcam-.CD35-.SLN", "LEC.SLN", "BEC.SLN", "IAP.SLN")

rna <- data.matrix(data.table::fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepers <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))

atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))

rna <- rna[,colnames(rna) %ni% blackout]
atac <- atac[keepers,colnames(atac) %ni% blackout]

rnaCPM <- sweep(rna, 2, colSums(rna), FUN="/") * 1000000
atacCPM <- sweep(atac, 2, colSums(atac), FUN="/") * 1000000

rnapc <- irlba(log2(rnaCPM+1), nv = 2)
atacpc <- irlba(log2(atacCPM+1), nv = 2)

rnadf <- data.frame(rnapc$v, names = colnames(rna)[colnames(rna) %ni% blackout])
atacdf <- data.frame(atacpc$v, names = colnames(atac)[colnames(atac) %ni% blackout])

ggplot(rnadf, aes(x = X1, y = X2, label = names)) + geom_text() + ggtitle("RNA")
ggplot(atacdf, aes(x = X1, y = X2, label = names)) + geom_text()+ ggtitle("ATAC")

