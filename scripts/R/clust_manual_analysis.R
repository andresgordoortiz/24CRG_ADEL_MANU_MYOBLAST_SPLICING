library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data (adjust file paths as needed)
clust_exons <- read.delim("notebooks/final/results/clust_results_tao_exons/Clusters_Objects.tsv",
                          header = TRUE, check.names = FALSE)[-1,]
clust_exons_matrix <- read.delim("notebooks/final/results/clust_results_tao_exons/Processed_Data/clust_events_splicing_tao_onlyexons.txt_processed.tsv",
                                 header = TRUE, check.names = FALSE)[-1,]

exons_longer <- pivot_longer(clust_exons_matrix, cols = -Genes, names_to = "time", values_to = "value")
exons_longer$significant <- exons_longer$Genes %in% differential_tao_exons_4d$EVENT

# Define your clusters; these should match the column names in clust_exons.
clusters <- c("C0 (672 genes)", "C1 (215 genes)", "C2 (250 genes)", "C3 (633 genes)")


# Compute percentage of significant events for each cluster
sig_data <- data.frame(Cluster = clusters, Percentage = NA)

for(i in seq_along(clusters)) {
  cl <- clusters[i]
  gene_list <- clust_exons[[cl]]
  
  percent_sig <- exons_longer %>% 
    filter(Genes %in% gene_list) %>%
    group_by(Genes) %>%
    summarise(isSig = any(significant)) %>%  # Check if any time point is significant
    summarise(pct = mean(isSig) * 100) %>%
    pull(pct)
  
  sig_data$Percentage[i] <- percent_sig
}

# Clean cluster names for the plot
sig_data$Cluster <- paste0("Cluster ", seq_along(clusters))

# Create bar plot
percent_plot <- ggplot(sig_data, aes(x = Cluster, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_brewer(palette = "Set2") +  # Professional color palette
  labs(title = "Percentage of Significant Exons per Cluster (4 Days vs 0 Days)",
       y = "Significant Events (%)",
       x = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

# Display the plot
print(percent_plot)


# Initialize an empty list to store the plots.
plots <- list()

# Loop over each cluster with an index to create new titles.
for(i in seq_along(clusters)){
  cl <- clusters[i]
  
  # Extract the list of genes for the current cluster.
  gene_list <- clust_exons[[cl]]
  
 
  # Extract the number of events from the original cluster string.
  event_count <- sub(".*\\((\\d+) genes\\).*", "\\1", cl)
  
  # Create a new title in the desired format: "Cluster i (event_count Events)"
  new_title <- paste0("Cluster ", i, " (", event_count, " Events)")
  
  # Create the plot for the current cluster.
  p <- ggplot(filter(exons_longer, Genes %in% gene_list), 
              aes(x = time, y = value, group = Genes, color = significant)) +
    geom_line(aes(alpha = significant)) + 
    scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "skyblue")) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
    labs(
      title = new_title,
      x = "Days of Differentiation",
      y = "Normalised PSI"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10),
      legend.position="none"
    )
  
  # Save the plot in the list.
  plots[[i]] <- p
}

# Combine the plots into a composite plot and add the big general title.
composite_plot_exons <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "Exon Splicing Patterns")

print(composite_plot_exons)



# Introns
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data (adjust file paths as needed)
clust_introns <- read.delim("notebooks/final/results/clust_results_tao_introns/Clusters_Objects.tsv",
                            header = TRUE, check.names = FALSE)[-1,]
clust_introns_matrix <- read.delim("notebooks/final/results/clust_results_tao_introns/Processed_Data/clust_events_splicing_tao_onlyintrons.txt_processed.tsv",
                                   header = TRUE, check.names = FALSE)[-1,]

introns_longer <- pivot_longer(clust_introns_matrix, cols = -Genes, names_to = "time", values_to = "value")
introns_longer$significant <- introns_longer$Genes %in% differential_tao_introns_4d$EVENT

# Define your clusters; these should match the column names in clust_introns.
clusters <- c("C0 (280 genes)", "C1 (266 genes)")

sig_data <- data.frame(Cluster = clusters, Percentage = NA)

for(i in seq_along(clusters)) {
  cl <- clusters[i]
  gene_list <- clust_introns[[cl]]
  
  percent_sig <- introns_longer %>% 
    filter(Genes %in% gene_list) %>%
    group_by(Genes) %>%
    summarise(isSig = any(significant)) %>%  # Check if any time point is significant
    summarise(pct = mean(isSig) * 100) %>%
    pull(pct)
  
  sig_data$Percentage[i] <- percent_sig
}

# Clean cluster names for the plot
sig_data$Cluster <- paste0("Cluster ", seq_along(clusters))

# Create bar plot
percent_plot_introns <- ggplot(sig_data, aes(x = Cluster, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_brewer(palette = "Set2") +  # Professional color palette
  labs(title = "Percentage of Significant Introns per Cluster (4 Days vs 0 Days)",
       y = "Significant Events (%)",
       x = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

# Display the plot
print(percent_plot_introns)

# Initialize an empty list to store the plots.
plots <- list()

# Loop over each cluster with an index to create new titles.
for(i in seq_along(clusters)){
  cl <- clusters[i]
  
  # Extract the list of genes for the current cluster.
  gene_list <- clust_introns[[cl]]
  
  # Extract the number of events from the original cluster string.
  event_count <- sub(".*\\((\\d+) genes\\).*", "\\1", cl)
  
  # Create a new title in the desired format: "Cluster i (event_count Events)"
  new_title <- paste0("Cluster ", i, " (", event_count, " Events)")
  
  # Create the plot for the current cluster.
  p <- ggplot(filter(introns_longer, Genes %in% gene_list), 
              aes(x = time, y = value, group = Genes, color = significant)) +
    geom_line(aes(alpha = significant)) + 
    scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "skyblue")) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
    labs(
      title = new_title,
      x = "Days of Differentiation",
      y = "Normalised PSI"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10),
      legend.position="none"
    )
  
  # Save the plot in the list.
  plots[[i]] <- p
}

# Combine the plots into a composite plot and add the big general title.
composite_plot_introns <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "Intron Splicing Patterns")

print(composite_plot_introns)

# Load required libraries
library(UpSetR)
library(tidyverse)

# Exons

# Define a function to extract significant events
extract_significant_events <- function(data) {
  filter(data, FDR <= 0.05 & abs(deltapsi) >= 0.1)$EVENT
}

# Create a named list of significant events per comparison
event_list <- rev(list(
  myo_2vs0 = extract_significant_events(tao_pdiff_exons_2d),
  myo_4vs0 = extract_significant_events(tao_pdiff_exons_4d),
  myo_7vs0 = extract_significant_events(tao_pdiff_exons)
))

# Convert to a binary presence/absence matrix
binary_matrix <- fromList(event_list)

# Generate an UpSet plot with improved aesthetics
UpSetR::upset(
  binary_matrix,
  sets = names(event_list),
  sets.bar.color = "steelblue",
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  keep.order = TRUE,
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)  # Adjust font sizes (default is 1)
)

# Introns

# Define a function to extract significant events
extract_significant_events <- function(data) {
  filter(data, FDR <= 0.05 & abs(deltapsi) >= 0.1)$EVENT
}

# Create a named list of significant events per comparison
event_list <- rev(list(
  myo_2vs0 = extract_significant_events(tao_pdiff_introns_2d),
  myo_4vs0 = extract_significant_events(tao_pdiff_introns_4d),
  myo_7vs0 = extract_significant_events(tao_pdiff_introns)
))

# Convert to a binary presence/absence matrix
binary_matrix <- fromList(event_list)

# Generate an UpSet plot with improved aesthetics
UpSetR::upset(
  binary_matrix,
  sets = names(event_list),
  sets.bar.color = "steelblue",
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  keep.order = TRUE,
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)  # Adjust font sizes (default is 1)
)
