library(irlba)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(cowplot)
"%ni%" <- Negate("%in%")

blackout <- c("Ep.MEChi.Th", "FRC.CD140a+.Madcam-.CD35-.SLN", "LEC.SLN", "BEC.SLN", "IAP.SLN")

rna <- data.matrix(data.table::fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))

rnapc <- irlba(log2(rna[,colnames(rna) %ni% blackout]+1), nv = 2)
atacpc <- irlba(log2(atac[,colnames(atac) %ni% blackout]+1), nv = 2)

rnadf <- data.frame(rnapc$v, names = colnames(rna)[colnames(rna) %ni% blackout])
atacdf <- data.frame(atacpc$v, names = colnames(atac)[colnames(atac) %ni% blackout])

ggplot(rnadf, aes(x = X1, y = X2, label = names)) + geom_text()
ggplot(atacdf, aes(x = X1, y = X2, label = names)) + geom_text()

