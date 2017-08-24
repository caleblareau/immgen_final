library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
"%ni%" <- Negate("%in%")

blackout <- c("Ep.MEChi.Th", "FRC.CD140a+.Madcam-.CD35-.SLN", "LEC.SLN", "BEC.SLN", "IAP.SLN")

# Import data
P <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_promoters.txt", sep = ","))
D <- data.matrix(read.table("../fromJDB/18-Aug-2017corr_distal-matchedP.txt", sep = ","))
sampleNames <- read.table("../fromJDB/18-Aug-2017corr.names.txt")[,1]
keep <- which(sampleNames %ni% blackout)
sampleNames <- sampleNames[keep]
P <- P[keep,keep]
D <- D[keep,keep]

# Import data
atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
atac <- atac[,sampleNames]

keepers <- data.table::fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepers <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))
peakAnno <- data.table::fread("../data/peakAnnoVector.txt")[[1]]
atac <- atac[1:dim(atac)[1] %in% keepers,]

# Variance component estimation
vals <- sapply(1:dim(atac)[1], function(i){
  Y <- log2(atac[i,] + 1)
  mod <- lmm.aireml(Y = scale(Y), K = list(P, D), verbose = FALSE)
  round(c(mod$sigma2, mod$tau),3)
})

# Make plot-ready data frame
vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
names(vdf) <- c("Unexplained", "Promoter", "Distal")
vdf$featureName <- peakAnno[ 1:length(peakAnno) %in% keepers]
vdf$peakID <-  data.table::fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]

write.table(vdf, file = "../output/peakVarianceComponents.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote= FALSE)

### Distal ####

vdf0 <- vdf[vdf$featureName == "outside",c(1,2,3)]
vdf0 <- vdf0[order(vdf0$Distal, decreasing = TRUE), ]
vdf2 <- rbind(vdf0[vdf0$Distal > 1, ],
              (vdf0[vdf0$Distal < 1, ])[(order(vdf0[vdf0$Distal < 1,"Promoter" ], decreasing = TRUE)),]
        )

vdf2 <- vdf2[sort(sample(1:dim(vdf2)[1],sum(vdf$featureName == "TSS"))),]
vdf2$rank <- 1:dim(vdf2)[1]
ldf <- melt(vdf2[,c("Unexplained", "Promoter", "Distal", "rank")], id.var=c("rank"))
p1 <- ggplot(ldf, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + pretty_plot() +
  labs(fill = "Component", x = "Distal Peaks", y = "% Variance Explained") + 
  scale_fill_manual(values = c("red", "green3", "dodgerblue")) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom") 

ggsave(p1, file = "../figures/23aug_varianceComponents_DISTAL.pdf", width = 10, height = 5)

vdf0 <- vdf[vdf$featureName == "TSS",]
vdf0 <- vdf0[order(vdf0$Promoter, decreasing = TRUE), ]
vdf2 <- rbind(vdf0[vdf0$Promoter > 1, ],
              (vdf0[vdf0$Promoter < 1, ])[(order(vdf0[vdf0$Promoter < 1,"Distal" ], decreasing = TRUE)),]
        )

vdf2$rank <- 1:dim(vdf2)[1]
ldf <- melt(vdf2[,c("Unexplained", "Promoter", "Distal", "rank")], id.var=c("rank"))
ldf$variable <- factor(as.character(ldf$variable), levels = c("Unexplained", "Distal", "Promoter"))
p1 <- ggplot(ldf, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + pretty_plot() +
  labs(fill = "Component", x = "Distal Peaks", y = "% Variance Explained") + 
  scale_fill_manual(values = c("red", "dodgerblue", "green3")) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom") 

ggsave(p1, file = "../figures/23aug_varianceComponents_PROMOTER.pdf", width = 10, height = 5)


