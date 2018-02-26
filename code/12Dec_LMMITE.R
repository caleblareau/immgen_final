library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(lmmlite)

# Import data
P <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_promoters.txt", sep = ","))
D <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_distal-matchedP.txt", sep = ","))
sampleNames <- read.table("../fromJDB/18-Aug-2017corr.names.txt")[,1]

rna <- data.matrix(data.table::fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
rna <- rna[,sampleNames]

geneNames <- read.table("../data/rnaRowLabels.txt", stringsAsFactors = FALSE)[,1]

SM <- read.table("../data/expressed_genes_SM.txt", stringsAsFactors = FALSE)[,1]

rna <- rna[geneNames %in% SM,]
geneNames <- geneNames[geneNames %in% SM]

# Variance component estimation-- compare GASTON to lmmlite
vals <- sapply(1:100, function(i){
  Y <- log2(rna[i,])
  mod <- lmm.aireml(Y = scale(Y), K = list(D), verbose = FALSE)
  round(c(mod$sigma2, mod$tau),3)
})
vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
names(vdf) <- c("Unexplained", "Distal")
vdf$gene <- geneNames[1:100]

vals2 <- sapply(1:100, function(i){
  Y <- log2(rna[i,])
  e <- eigen_rotation(D, Y, rep(1, length(Y)))
  fitLMM(e$Kva, e$y, e$X)$hsq
})
vdf$lmmite <- vals2

##
# okay so this looks good
##

vals <- sapply(1:dim(rna)[1], function(i){
  Y <- log2(rna[i,])
  mod <- lmm.aireml(Y = scale(Y), K = list(D,P), verbose = FALSE)
  round(c(mod$sigma2, mod$tau),3)
})
vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
names(vdf) <- c("Unexplained", "Distal", "Promoter")
vdf$gene <- geneNames
