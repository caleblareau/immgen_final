library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(varistran)
library(lmmlite)
load("usefulData/everythingForVC_CL.rda")

P <- cor(promoters)
D <- cor(all_distal)

y <- varistran::vst(rna, cpm = TRUE)
rna <- y

permuteSeveral <- function(i){
  set.seed(i)
  rna2 <- rna[,sample(1:80)]

  # Variance component estimation
  vals <- sapply(1:dim(rna)[1], function(i){
    Y <- scale(rna2[i,])
    mod <- lmm.aireml(Y = Y, K = list(P, D), verbose = FALSE)
    round(c(mod$sigma2, mod$tau),3)
  })
  
  # Add in the p-value from lmmlite
  #pvals <- sapply(1:dim(rna)[1], function(i){
  #  Y <- (rna[i,])
  #  e <- eigen_rotation(D, Y, matrix(rep(1, length(Y), ncol = 1)))
  #  out <- fitLMM(Kva = e$Kva, y = e$y, X = e$X, compute_se = TRUE, use_cpp=FALSE)
  #})
  
  
  
  # Make plot-ready data frame
  vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
  names(vdf) <- c("Unexplained", "Promoter", "Distal")
  vdf$gene <- rownames(rna)
  
  write.table(vdf, file = paste0("permute_output/vco-p",as.character(i),".txt"),
                                 quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
lapply(1:100, permuteSeveral)



