library(ggplot2)
library(BuenColors)

distalDF <- head(data.frame(
  PATHWAY = c("immune system process", "regulation of immune system process", "positive regulation of immune system process", 
              "signal transduction", "regulation of cell adhesion", "regulation of cell activation", "cell activation", 
              "positive regulation of biological process", "response to stimulus", "immune response"),
  FDR = c(3E-27, 1E-24, 2E-21, 1.88E-20, 7.8E-20, 5E-19, 1E-18, 2E-18, 2E-17, 2E-17)
),5)

promoterDF <- data.frame(
  PATHWAY = c("peptide biosynthetic process", "amide biosynthetic process", "translation", "organonitrogen compound biosystenic process", "DNA replication"),
  FDR = c(8E-8, 7E-8, 2E-7, 2E-7, 6E-6)
)

unexplainedDF <- data.frame(
  PATHWAY = c("system process", "sensory perception", "sensotry perception of chemical stimulus", "nervous system process", "G-protein coupled receptor signaling pathway"),
  FDR = c(2E-4, 8E-4, 3E-3, 2E-3, 2E-3)
)

plotdf <- rbind(distalDF, promoterDF, unexplainedDF)
plotdf$Component <- c(rep("Distal", dim(distalDF)[1]), rep("Promoter", dim(promoterDF)[1]),
                 rep("Unexplained", dim(unexplainedDF)[1]))
plotdf$logFDR <- -1*log10(plotdf$FDR)
plotdf$PATHWAY <- as.character(plotdf$PATHWAY)
plotdf$PATHWAY <- factor(plotdf$PATHWAY, levels = plotdf$PATHWAY )
  
pp <- ggplot(data=plotdf, aes(x=PATHWAY, y=logFDR, fill = Component)) +
  geom_bar(stat="identity", color = "black") +
  scale_fill_manual(values = c("dodgerblue", "green3", "red")) +
  coord_flip() + pretty_plot() + theme( legend.position="bottom") +
  labs(fill = "Variance component", y = "- log10 FDR", x = "Pathway") +
  ggtitle("Gene Ontology Term Enrichment")
ggsave(pp, file = "../figures/24august_VC_GO_enrich.pdf")
