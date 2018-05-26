library(data.table)
library(dplyr)

doImport <- function(i){
  file <- paste0("permute_output/vco-p",as.character(i),".txt")
  tab <- read.table(file, header = TRUE)
  tab$i <- i
  tab
}

lapply(1:16, doImport) %>% rbindlist() %>% as.data.frame() -> df

df %>% group_by(gene) %>% summarise(meanPermute_Unexplained = mean(Unexplained), 
                                    meanPermute_Promoter = mean(Promoter), 
                                    meanPermute_Distal = mean(Distal)) %>% as.data.frame() -> mean_permuted_df
                 
observed <- read.table("output/varianceComponentsOutput.txt", header = TRUE)
mdf <- merge(mean_permuted_df, observed)

mdf$density <- get_density(mdf$meanPermute_Unexplained, mdf$Distal)
ggplot(mdf, aes(meanPermute_Unexplained, Distal, color = density)) + geom_point() +
  labs(x = "Mean Permuted % Unexplained", y = "True Observed Distal % Explained") +
  ggtitle("Each dot is a gene") + pretty_plot() + L_border()

cor(mdf$meanPermute_Unexplained, mdf$Unexplained)
