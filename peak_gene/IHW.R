library(IHW)
library(ggplot2)

dt <- data.table::fread("zcat < Peak-to-gene_perPeak_Pearson_noStromal.txt.gz", header = FALSE, sep = "\t")
dt$distGroup <- groups_by_filter(abs(dt[["V3"]]), 6)
dt$p.value <- 10^(-1*dt[["V4"]])

ihwRes <- ihw(p.value ~ distGroup, data = dt, alpha = 0.1)

ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)

