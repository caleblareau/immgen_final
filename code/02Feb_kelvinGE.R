library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)

rna <- data.matrix(data.table::fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
geneNames <- read.table("../data/rnaRowLabels.txt", stringsAsFactors = FALSE)[,1]
rownames(rna) <- geneNames

log2cpm <- log2(sweep(rna, 2, colSums(rna), FUN="/") * 1000000 + 1)

qplot(rowMeans(log2cpm), bins = 100) + pretty_plot() 
