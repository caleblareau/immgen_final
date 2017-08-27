library(IHW)
library(ggplot2)
library(BuenColors)

dt <- data.table::fread("zcat < Peak-to-gene_perPeak_Pearson_noStromal.txt.gz", header = FALSE, sep = "\t")
dt$distGroup <-cut(dt[["V3"]],seq(-1000000,1000000,by = 250000))
dt <- dt[complete.cases(dt),]
dt$p.value <- 10^(-1*dt[["V4"]])

ihwRes <- ihw(p.value ~ distGroup, data = dt, alpha = 0.1)
ihwResdf <- as.data.frame(ihwRes)
ihwResdf$padjBH <- p.adjust(dt[["p.value"]], method = "BH")
ihwResdf$logFC <- log10(log10(ihwResdf$adj_pvalue)/log10(ihwResdf$padjBH))
ihwResdf$distance <- dt[["V3"]]

p1 <- ggplot(dt, aes(x=p.value)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ distGroup, nrow = 3) + pretty_plot()
ggsave(p1, file = "stratifiedPvalueHistogram.pdf")

p2 <- ggplot(ihwResdf, aes(x=group, y=logFC)) + 
  geom_boxplot()   + pretty_plot()   
ggsave(p2, file = "stratifiedBoxPlotFoldChange.pdf")

#ggplot(dt, aes(x = p.value, col = distGroup)) + stat_ecdf(geom = "step")+ pretty_plot()
