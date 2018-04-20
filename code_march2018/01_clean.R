library(data.table)
"%ni%" <- Negate("%in%")

rna <- data.matrix(data.table::fread("cat < usefulData/MF.fixed.RNAseq.pop.mean.20180213csv.csv"))[,-1]
genes <- data.frame(data.table::fread("cat < usefulData/MF.fixed.RNAseq.pop.mean.20180213csv.csv"))[,1]

stromal <- c("FRC.CD140a+.Madcam-.CD35-.SLN", "LEC.SLN", "IAP.SLN", "BEC.SLN", "Ep.MEChi.Th")
rna <- rna[,colnames(rna) %ni% c(stromal, "MF.ICAM+480hi.PC")]
rownames(rna) <- genes


# Import ATAC and filter peaks
atac <- data.matrix(data.table::fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))
keepers <- data.table::fread("../data/immgen_good_peaks.txt", header = FALSE)[[1]]
keepers <-  as.numeric(gsub("ImmGenATAC1219.peak_", "", keepers))
peakAnno <- data.table::fread("../data/peakAnnoVector.txt")[[1]]
atac <- atac[1:dim(atac)[1] %in% keepers,]
peakAnno <- peakAnno[1:length(peakAnno) %in% keepers]

# Filter out samples such that they are all in the RNA
colnames(rna) [colnames(rna) %ni% colnames(atac)]
atac <- atac[,colnames(rna)]

# Filter Sara's genes
Exg <- read.table("usefulData/expressed_genes_SM.txt", stringsAsFactors = FALSE)[,1]
Exg <- Exg[Exg %in% rownames(rna)]
rna <- rna[Exg,]
#rna <- log2(rna + 1)

promoters <- atac[peakAnno == "TSS",]
all_distal <- atac[peakAnno == "outside",]

save(promoters, all_distal, rna, file = "usefulData/everythingForVC_CL.rda")


