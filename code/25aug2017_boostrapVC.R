library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)

# Import data
P <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_promoters.txt", sep = ","))
D <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_distal-matchedP.txt", sep = ","))
sampleNames <- read.table("../fromJDB/18-Aug-2017corr.names.txt")[,1]

rna <- data.matrix(data.table::fread("zcat < ../data/17aug_commonRNA_normalized.txt.gz"))
rna <- rna[,sampleNames]
bsout <- lapply(1:100, function(k){
  randVec <- sample(1:dim(P)[1], dim(P)[1], replace = TRUE)
  # Variance component estimation
  vals <- sapply(1:dim(rna)[1], function(i){
    Y <- log2(rna[i,randVec])
    mod <- lmm.aireml(Y = scale(Y), K = list(P[randVec,randVec], D[randVec,randVec]), verbose = FALSE)
    round(c(mod$sigma2, mod$tau),3)
  })
  
  # Make plot-ready data frame
  vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
  names(vdf) <- c("Unexplained", "Promoter", "Distal")
  vdf
})
saveRDS(bsout, file = "bootstrapRaw.rds")
