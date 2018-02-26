library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(irlba)

# Import data
P <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_promoters.txt", sep = ","))
D <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_distal-matchedP.txt", sep = ","))
sampleNames <- read.table("../fromJDB/18-Aug-2017corr.names.txt")[,1]

rna <- data.matrix(data.table::fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
rna <- rna[,sampleNames]
atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
atac <- atac[,sampleNames]
geneNames <- read.table("../data/rnaRowLabels.txt", stringsAsFactors = FALSE)[,1]

atacpcs <- prcomp_irlba(x = atac, n = 20)
pcdf <- data.frame(atacpcs$rotation)

r2 <- sapply(1:length(geneNames), function(i){
  print(i)
  mod <- lm(rna[i,] ~ pcdf$PC1 +  pcdf$PC2 +  pcdf$PC3 + pcdf$PC4 + pcdf$PC5 + pcdf$PC6 +
              pcdf$PC7 + pcdf$PC8 + pcdf$PC9 + pcdf$PC10 + pcdf$PC11 + pcdf$PC12 + pcdf$PC13 + pcdf$PC14 +  pcdf$PC15 + 
              pcdf$PC16 + pcdf$PC17 + pcdf$PC18 + pcdf$PC19 + pcdf$PC20)
  return(c(summary(mod)$r.squared, summary(mod)$adj.r.squared))
})

r2df <- data.frame(geneNames, t(r2))

vcdf <- data.frame(data.table::fread("../output/RNA_varianceComponents.txt", header = TRUE))

mdf <- merge(vcdf, r2df, by.x = "gene", by.y = "geneNames")
pdf("../figures/PCregression.pdf", height = 5, width = 5)
smoothScatter(mdf$Unexplained, (1- mdf$X2)*100, xlab = "% Variance Unexplained in Variance Component Model", ylab = "% Variance Unexplained in PC Regression")
dev.off()
cor(mdf$Unexplained, 1 - mdf$X2)


