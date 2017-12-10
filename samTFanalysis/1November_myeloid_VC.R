library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(readxl)
library(igraph)
library(data.table)
"%ni%" <- Negate("%in%")

#############
# ATAC zscores
##############

dev <- read.csv("../chromVAR/deviationsData.csv", header = FALSE)
colnames(dev) <- read.csv("../chromVAR/deviationsCols.tsv", header = FALSE)[,1]
row.names(dev) <- read.csv("../chromVAR/deviationsRows.tsv", header = FALSE)[,1]

# Import Sam data
mat <- (readRDS("data/correlationDistanceLog2_myeloid_noProgenitors.rds"))
labs <- attributes(mat)$Labels

edges <- data.frame(read_excel("../data/PopulationTree.xls")[,c(1,2)])
colnames(edges) <- c("source", "target")
ig <- graph.data.frame(edges)
s.paths <- shortest.paths(ig, algorithm = "dijkstra")

all <- intersect(colnames(dev), colnames(s.paths)) %>% intersect(labs)
length(all)
labs[labs %ni% all]

colnames(s.paths)[!(colnames(s.paths) %in% all)]
colnames(dev)[!(colnames(dev) %in% all)]

# Update to those that intersect
related <- 0.5^s.paths[all, all]
devmodel <- data.matrix(setcolorder(dev[,colnames(dev) %in% all],  all))

dim(related)
pearson <- 1- as.matrix(mat)[all,all] 

vals <- sapply(1:dim(devmodel)[1], function(i){
  Y <- (devmodel[i,])
  mod <- lmm.aireml(Y = scale(Y), K = list(pearson), verbose = FALSE)
  round(c(mod$sigma2, mod$tau),3)
})

vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
names(vdf) <- c("Unexplained", "Phylogeny")
vdf$TF <- row.names(dev)
head(vdf)

vdf <- vdf[order(vdf$Phylogeny, decreasing = TRUE), ]
vdf2 <- vdf
vdf2$rank <- 1:dim(vdf2)[1]

ldf <- melt(vdf2, id.var=c("TF", "rank"))
ggplot(ldf, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + pretty_plot() +
  labs(fill = "Component", x = "TF", y = "% Variance Explained") + 
  scale_fill_manual(values = c("purple4", "orange")) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom") 

write.table(vdf2, file = "myeloid_lineageDefined_1Novem.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(vdf2, file = "../output/25september_tissuePhylogeny_VC.tsv", quote = FALSE, sep = "\t",
            col.names = TRUE, row.names = FALSE)

#########
# Peaks
#########

atacRaw <- data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz")
keepers <- fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepers2 <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))
atacRaw <- atacRaw[keepers2,]

edges <- data.frame(read_excel("../data/PopulationTree.xls")[,c(1,2)])
colnames(edges) <- c("source", "target")
ig <- graph.data.frame(edges)
s.paths <- shortest.paths(ig, algorithm = "dijkstra")

all <- intersect(colnames(atacRaw), colnames(s.paths))
length(all)

colnames(s.paths)[!(colnames(s.paths) %in% all)]
colnames(rna)[!(colnames(rna) %in% all)]

# Update to those that intersect
related <- 0.5^s.paths[all, all]
atacmodel <- data.matrix(setcolorder(atacRaw[,colnames(atacRaw) %in% all, with=FALSE],  all))

dim(related)
dim(atacmodel)

# Make tissue covariance matrix
tissues <- sapply(strsplit(colnames(atacmodel), split = "[.]"), function(i) i[length(i)])
tissues[tissues == "18hr"] <- "Sp"
tissues[tissues == "3hr"] <- "Sp"
tissues[tissues == "3d"] <- "Sp"
table(tissues)
Tissue <- sapply(tissues, function(t) as.numeric(t == tissues))

vals <- sapply(1:dim(atacmodel)[1], function(i){
  Y <- log2(atacmodel[i,])
  mod <- lmm.aireml(Y = scale(Y), K = list(Tissue, related), verbose = FALSE)
  round(c(mod$sigma2, mod$tau),3)
})

vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
names(vdf) <- c("Neither", "Tissue", "Phylogeny")
head(vdf)

vdf <- vdf[order(vdf$Phylogeny, decreasing = TRUE), ]
vdf2 <- rbind(vdf[vdf$Phylogeny > 1, ],
              (vdf[vdf$Phylogeny < 1, ])[(order(vdf[vdf$Phylogeny < 1,"Tissue" ], decreasing = TRUE)),]
        )
vdf2$rank <- 1:dim(vdf2)[1]

ldf <- melt(vdf2, id.var=c("rank"))
p1 <-ggplot(ldf, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + pretty_plot() +
  labs(fill = "Component", x = "Genes", y = "% Variance Explained") + 
  scale_fill_manual(values = c("purple4", "orange", "yellow")) + ggtitle("Peaks in Celltypes / Phylogeny") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom") 
ggsave(p1, file = "../figures/26September_varianceComponents2Peaks.pdf")

write.table(vdf2, file = "../output/25september_tissuePhylogeny_VC.tsv", quote = FALSE, sep = "\t",
            col.names = TRUE, row.names = FALSE)

