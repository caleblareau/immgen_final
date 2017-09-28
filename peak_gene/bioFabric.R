library(igraph)
library(ggnet)

dt <- data.table::fread("zcat < Peak-to-gene_perPeak_Pearson_noStromal.txt.gz", header = FALSE, sep = "\t")
ig <- graph_from_data_frame(dt[,c(1,2)])
adj <- as_adj(ig, type = c("both"))
adj <- network(adj, directed = FALSE)

