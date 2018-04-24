library(gaston)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(varistran)
library(mgatk)
library(matrixStats)

load("usefulData/everythingForVC_CL.rda")

#y <- varistran::vst(rna, cpm = TRUE)
#rna <- y
rna <- log2(rna+1)
g <- giniRows(rna)
CV <- sqrt(rowVars(rna))/rowMeans(rna)
Exp <- rowMeans(rna)

vco <- data.frame(data.table::fread("output/varianceComponentsOutput.txt"))
thres <- 99
dg <- as.character(vco$gene)[vco$Distal > thres] 
pg <- as.character(vco$gene)[vco$Promoter > thres] 
ug <- as.character(vco$gene)[vco$Unexplained > thres] 

write.table(data.frame(pg), quote = FALSE, sep  = "\t", row.names = FALSE, col.names = FALSE)

gdf <- data.frame(
  Gini = c(g[rownames(rna) %in% dg], g[rownames(rna) %in% pg], g[rownames(rna) %in% ug]),
  CV = c(CV[rownames(rna) %in% dg], CV[rownames(rna) %in% pg], CV[rownames(rna) %in% ug]),
  Exp = c(Exp[rownames(rna) %in% dg], Exp[rownames(rna) %in% pg], Exp[rownames(rna) %in% ug]),
  Type = c(rep("Distal", sum(rownames(rna) %in% dg)), rep("Promoter", sum(rownames(rna) %in% pg)), rep("Unexplained",sum(rownames(rna) %in% ug)))
)

p1 <- ggplot(gdf, aes(x = Type, y = Gini, color = Type)) + geom_boxplot() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "", y = "Gini",  color = "") +
  scale_color_manual(values = c("dodgerblue3", "green3", "red"))

p2 <- ggplot(gdf, aes(x = Type, y = CV, color = Type)) + geom_boxplot() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "", y = "Coefficient of Variation",  color = "") +
  scale_color_manual(values = c("dodgerblue3", "green3", "red"))

p3 <- ggplot(gdf, aes(x = Type, y = Exp, color = Type)) + geom_boxplot() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "", y = "Average log2 Expression",  color = "") +
  scale_color_manual(values = c("dodgerblue3", "green3", "red"))


cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, nrow = 1), filename = "figures/FigureS2.pdf", width = 12, height = 2.5)
