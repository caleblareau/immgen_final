
library(data.table)

genedf <- fread("../data/mm10.refGenes.2016.1018.csv")
vcdf <- fread("../output/19aug_varianceComponents.txt")
load("../mouseGSEA/immgenMSIGDBlists.rda")

enrich<- function(geneVec, pathwayList, total = length(promoterGenes) + length(distalGenes) + length(unexplainedGenes)){
  hg_pval_hits_genes <- function(pathwayVec, genesIn, total){
    hits <- sum(genesIn %in% pathwayVec)
    genes <- paste(genesIn[genesIn %in% pathwayVec], collapse = ",")
    return(c(as.character(1-phyper(hits - 1, length(pathwayVec), total - length(pathwayVec), length(genesIn))), 
                      as.character(hits), genes))
  }
  
  df <- data.frame(t(as.data.frame(lapply(pathwayList, hg_pval_hits_genes, genesIn = geneVec, total = total))))
  df$Pathway <- rownames(df)
  rownames(df) <- NULL
  colnames(df) <- c("Pvalue", "Hits", "Genes", "Pathway")
  df$FDR <- p.adjust(df$Pvalue, method = "BH")
  return(df[order(df$Pvalue),c(4,5,1,2,3)])
}

promoterGenes <- vcdf[vcdf$Promoter > 99.99,][["gene"]]
distalGenes <- vcdf[vcdf$Distal > 99.99,][["gene"]]
unexplainedGenes <- vcdf[vcdf$Unexplained > 99.99,][["gene"]]

length(promoterGenes)
length(distalGenes)
length(unexplainedGenes)

head(enrich(promoterGenes, mm.H)[,1:4])
head(enrich(distalGenes, mm.H)[,1:4])
head(enrich(unexplainedGenes, mm.H)[,1:4])

head(enrich(promoterGenes, mm.c2)[,1:4])
head(enrich(distalGenes, mm.c2)[,1:4])
head(enrich(unexplainedGenes, mm.c2)[,1:4])

head(enrich(promoterGenes, mm.c3)[,1:4])
head(enrich(distalGenes, mm.c3)[,1:4])
head(enrich(unexplainedGenes, mm.c3)[,1:4])

head(enrich(promoterGenes, mm.c4)[,1:4])
head(enrich(distalGenes, mm.c4)[,1:4])
head(enrich(unexplainedGenes, mm.c4)[,1:4])

head(enrich(promoterGenes, mm.c5)[,1:4])
head(enrich(distalGenes, mm.c5)[,1:4])
head(enrich(unexplainedGenes, mm.c5)[,1:4])

head(enrich(promoterGenes, mm.c6)[,1:4])
head(enrich(distalGenes, mm.c6)[,1:4])
head(enrich(unexplainedGenes, mm.c6)[,1:4])

head(enrich(promoterGenes, mm.c7)[,1:4])
head(enrich(distalGenes, mm.c7)[,1:4])
head(enrich(unexplainedGenes, mm.c7)[,1:4])
