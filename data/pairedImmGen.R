library(data.table)

atac <- data.matrix(fread("zcat < ATAC.density.normalized.csv.gz", header = TRUE)[,-1])
rna <- data.matrix(fread("zcat < RNAseq.pop.mean.csv.gz", header = TRUE)[,-1])
common <- intersect(colnames(rna),colnames(atac))

c.atac <- round(atac[,common],2)
c.rna <- round(rna[,common],2)
all(colnames(c.rna) == colnames(c.atac))

write.table(c.atac, file = "17aug_commonATAC_normalized.txt", sep = "\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE)

write.table(c.rna, file = "17aug_commonRNA_normalized.txt", sep = "\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE)