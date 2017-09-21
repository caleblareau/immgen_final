library(data.table)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifmatchr)
library(chromVARmotifs)
library(GenomicRanges)
library(BuenColors)

# Import human data
hemecounts <- data.matrix(fread("zcat < ../humandata/GSE74912_ATACseq_All_Counts.txt.gz"))
h <- hemecounts[,4:80]
 cc <- colnames(h)
 shortnames <- sapply(1:length(cc), function(s){
   if(s <= 23 | (s <= 74 & s >= 62 )){ strsplit(cc[s], split = "-")[[1]][2]
   } else { strsplit(cc[s], split = "-")[[1]][3]}
 })
noa <- gsub("A", "", shortnames)
nob <- gsub("B", "", noa)
noc <- gsub("C", "", nob)
#
n <- as.character(seq(1:15))
cell <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Mono", "missing", "CD4", "CD8", "NK", "missing", "Bcell", "CLP", "Ery")
names(cell) <- n
namess <- paste0(unname(cell[noc]))
colnames(h) <- namess

cellTypeDatMatHuman <- sapply(unique(namess), function(namee){
  rowSums(h[,c(which(namee== namess))])
})
remove(hemecounts)

humanPeaks <- fread("../humandata/peaks.bed")
humanPeaks$idx <- paste0("human", as.character(1:dim(humanPeaks)[1]))
human_g <-  makeGRangesFromDataFrame(humanPeaks, seqnames = "Chr", start.field = "Start",
                                     end.field = "End", keep.extra.columns = TRUE)

remove(h)
counts <- SummarizedExperiment(assays = list(counts = cellTypeDatMatHuman),
                               rowData = human_g, 
                               colData = DataFrame(names = colnames(cellTypeDatMatHuman)))
counts <- filterPeaks(counts)

counts <- addGCBias(counts,  genome = BSgenome.Hsapiens.UCSC.hg19)

makeHitMat <- function(file){
  gwas_g <- makeGRangesFromDataFrame(fread(file), seqnames = "V1", start.field = "V2", end.field = "V3")
  as.numeric(1:length(counts@rowRanges) %in% queryHits(findOverlaps(counts@rowRanges, gwas_g)))
}

files <- list.files("../gwasRaw", full.names = TRUE)
names <- gsub("../gwasRaw/", "", files) %>% gsub(pattern = "_pruned_rsq_0.8.bed", replacement = "")

hitMat250 <- sapply(files, makeHitMat)
colnames(hitMat250) <- (names)
hitMat250 <- Matrix::Matrix(hitMat250[,colSums(hitMat250) >= 50], sparse = TRUE)

dev <- computeDeviations(object = counts, annotations = hitMat250)
gwasz <- reshape2::melt(assays(dev)[["z"]])
sorted_gwasz <- gwasz[order(gwasz$value, decreasing = TRUE),]

