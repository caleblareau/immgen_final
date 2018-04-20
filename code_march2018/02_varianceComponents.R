library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(varistran)
library(lmmlite)

load("usefulData/everythingForVC_CL.rda")

P <- cor(promoters)
D <- cor(all_distal)
A <- cor(atac)

y <- varistran::vst(rna, cpm = TRUE)
rna <- y

# Variance component estimation
vals <- sapply(1:dim(rna)[1], function(i){
  Y <- scale(rna[i,])
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
vdf <- vdf[order(vdf$Distal, decreasing = TRUE), ]
vdf2 <- rbind(vdf[vdf$Distal > 1, ],
              (vdf[vdf$Distal < 1, ])[(order(vdf[vdf$Distal < 1,"Promoter" ], decreasing = TRUE)),]
        )
vdf2$rank <- 1:dim(vdf2)[1]

ldf <- melt(vdf2, id.var=c("gene", "rank"))
p1 <- ggplot(ldf, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + pretty_plot(fontsize = 6) +
  labs(fill = "Component", x = "Genes", y = "% Variance Explained") + 
  scale_fill_manual(values = c("red", "green3", "dodgerblue3")) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "none") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_x_continuous( expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 

cowplot::ggsave(p1, file = "figures/newVarianceComponent.pdf", height = 2, width = 4)

write.table(vdf2, file = "output/varianceComponentsOutput.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)




  
  