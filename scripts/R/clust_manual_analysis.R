library(tidyverse)
clust_exons<-read.delim("notebooks/final/results/clust_results_tao_exons/Clusters_Objects.tsv", header = TRUE, check.names = FALSE)[-1,]
clust_exons_matrix<-read.delim("notebooks/final/results/clust_results_tao_exons/Processed_Data/clust_events_splicing_tao_onlyexons.txt_processed.tsv", header = TRUE, check.names = FALSE)[-1,]

exons_longer<-pivot_longer(clust_exons_matrix, cols = -Genes, names_to = "time", values_to = "value")
exons_longer$significant<-exons_longer$Genes %in% differential_tao_exons_4d$EVENT


library(dplyr)
library(ggplot2)
library(patchwork)

# Define your clusters; these should match the column names in clust_exons.
clusters <- c("C0 (672 genes)", "C1 (215 genes)", "C2 (250 genes)", "C3 (633 genes)")

# Initialize an empty list to store the plots.
plots <- list()

# Loop over each cluster
for(cl in clusters){
  
  # Extract the list of genes for the current cluster.
  gene_list <- clust_exons[[cl]]
  
  # Calculate percentage of significant genes for the current cluster.
  percent_sig <- exons_longer %>% 
    filter(Genes %in% gene_list) %>%
    group_by(Genes) %>%
    summarise(isSig = any(significant)) %>%  # if any time point is significant for that gene
    summarise(pct = mean(isSig) * 100) %>%
    pull(pct)
  
  # Create the plot for the current cluster.
  p <- ggplot(filter(exons_longer, Genes %in% gene_list), 
              aes(x = time, y = value, group = Genes, color = significant)) +
    geom_line(aes(alpha = significant)) + 
    scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "skyblue")) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
    labs(
      title = paste("Gene Expression Trends for", cl),
      subtitle = paste0("Percentage of Significant Genes: ", round(percent_sig, 1), "%"),
      x = "Time",
      y = "Expression Value",
      caption = "Data source: exons_longer"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10)
    )
  
  # Save the plot in the list.
  plots[[cl]] <- p
}

# Combine the plots into a composite plot.
composite_plot <- wrap_plots(plots, ncol = 2)
print(composite_plot)



