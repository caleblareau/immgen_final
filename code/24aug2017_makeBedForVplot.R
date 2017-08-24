
# Make .bed files for the 

vcdf <- read.table("../output/19aug_varianceComponents.txt", header = TRUE)
gcdf <- read.table("../data/mm10.refGenes.2016.1018.csv", sep = ",", header = TRUE)

promoterGenes <- vcdf[vcdf$Promoter > 99.9,"gene"]
distalGenes <- vcdf[vcdf$Distal > 99.9,"gene"]
unexplainedGenes <- vcdf[vcdf$Unexplained > 99.9,"gene"]

write.table(gcdf[gcdf[,3] %in% promoterGenes,c(4,8,8,5)], file = "../output/vplot_bed/promotersVC.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(gcdf[gcdf[,3] %in% distalGenes,c(4,8,8,5)], file = "../output/vplot_bed/distalVC.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(gcdf[gcdf[,3] %in% unexplainedGenes,c(4,8,8,5)], file = "../output/vplot_bed/unexplainedVC.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

