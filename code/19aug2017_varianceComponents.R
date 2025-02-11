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
atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
atac <- atac[,sampleNames]
geneNames <- read.table("../data/rnaRowLabels.txt", stringsAsFactors = FALSE)[,1]

# Variance component estimation
vals <- sapply(1:dim(rna)[1], function(i){
  Y <- log2(rna[i,])
  mod <- lmm.aireml(Y = scale(Y), K = list(P, D), verbose = FALSE)
  round(c(mod$sigma2, mod$tau),3)
})


# Make plot-ready data frame
vdf <- data.frame(t(vals)/rowSums(t(vals))*100)
names(vdf) <- c("Unexplained", "Promoter", "Distal")
vdf$gene <- geneNames
vdf <- vdf[order(vdf$Distal, decreasing = TRUE), ]
vdf2 <- rbind(vdf[vdf$Distal > 1, ],
              (vdf[vdf$Distal < 1, ])[(order(vdf[vdf$Distal < 1,"Promoter" ], decreasing = TRUE)),]
        )
vdf2$rank <- 1:dim(vdf2)[1]


ldf <- melt(vdf2, id.var=c("gene", "rank"))
p1 <- ggplot(ldf, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + pretty_plot() +
  labs(fill = "Component", x = "Genes", y = "% Variance Explained") + 
  scale_fill_manual(values = c("red", "green3", "dodgerblue")) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom") 
ggsave(p1, file = "../figures/20aug_varianceComponents.pdf")

# Import the peak annotation; filter things that are not in 
# the RNA-seq data or duplicated TSSs
annoDF <- data.table::fread("../data/annotated.peaks.with.myID.txt")
TSSdf <- annoDF[annoDF[[6]] == "TSS",]
dTSSdf <- TSSdf[!duplicated(TSSdf$genename.of.closest.TSS), ]
TSSdfFinal <- dTSSdf[dTSSdf[["genename.of.closest.TSS"]] %in% vdf2$gene,]

# Do promoter - gene association
totalDFlist <- lapply(1:dim(TSSdfFinal)[1], function(i){
  atacidx <- as.numeric(strsplit(TSSdfFinal[i,][["peak.ID.number"]], split = "_")[[1]][2])
  gene <-  TSSdfFinal[i,][["genename.of.closest.TSS"]]
  vcidx <- which(vdf2$gene == gene)
  df <- data.frame(vdf2[vcidx,], peakID = TSSdfFinal[i,][["peak.ID.number"]])
  df$pearson <- cor(log2(rna[which(geneNames==gene), ]), log2(atac[atacidx,]))
  df
})

totalDF <- data.table::rbindlist(totalDFlist)

vcPeak <- data.frame(data.table::fread("../output/peakVarianceComponents.tsv"))
totalDF2 <- merge(totalDF, vcPeak[,c(2,5)], by = "peakID")

#write.table(totalDF, file = "../output/19aug_varianceComponentsWpearson.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

totalDF2$Bins <- cut(totalDF2$Promoter.x, breaks=seq(0,100,10))
ggplot(totalDF2, aes(x = Bins, y = pearson^2)) + 
  geom_boxplot()+ pretty_plot()



  
  