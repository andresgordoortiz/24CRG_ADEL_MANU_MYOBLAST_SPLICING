library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)  # for refined color palette

# Load data (adjust file paths as needed)
clust_exons <- read.delim("notebooks/final/results/clust_results_tao_exons/Clusters_Objects.tsv",
                          header = TRUE, check.names = FALSE)[-1,]
clust_exons_matrix <- read.delim("notebooks/final/results/clust_results_tao_exons/Processed_Data/clust_events_splicing_tao_onlyexons.txt_processed.tsv",
                                 header = TRUE, check.names = FALSE)[-1,]

# Convert to long format
exons_longer <- pivot_longer(clust_exons_matrix, cols = -Genes, names_to = "time", values_to = "value")

# Define your clusters; these names should match the column names in clust_exons.
clusters <- c("C0 (672 genes)", "C1 (215 genes)", "C2 (250 genes)", "C3 (633 genes)")

# Create a mapping from the original cluster names to new labels.
cluster_map <- data.frame(
  Original = clusters,
  New = sapply(seq_along(clusters), function(i) {
    event_count <- length(clust_exons[[clusters[i]]])
    paste0("Cluster ", i, " (", event_count, " Events)")
  }),
  stringsAsFactors = FALSE
)

# Create a list of differential comparisons
comparisons <- list(
  "4d" = differential_tao_exons_4d,
  "7d" = differential_tao_exons,
  "2d" = differential_tao_exons_2d,
  "4dvs2d" = differential_tao_exons_4d_2d
)

# Initialize an empty data frame for counts per comparison and cluster.
stacked_data <- data.frame(Comparison = character(),
                           Cluster = character(),
                           sig_count = numeric(),
                           stringsAsFactors = FALSE)

# Loop over comparisons and clusters.
for(comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  
  # Mark significance based on the current differential result.
  exons_longer$significant <- exons_longer$Genes %in% comp$EVENT
  
  for(cl in clusters) {
    gene_list <- clust_exons[[cl]]
    
    # Count unique genes in the cluster that are significant.
    count_sig <- exons_longer %>% 
      filter(Genes %in% gene_list) %>%
      group_by(Genes) %>%
      summarise(isSig = any(significant)) %>%
      summarise(count = sum(isSig)) %>%
      pull(count)
    
    # In case there are no significant events, set count to zero.
    if(length(count_sig) == 0) count_sig <- 0
    
    # Get the new cluster label.
    new_cluster <- cluster_map$New[which(cluster_map$Original == cl)]
    
    # Append the result for this comparison and cluster.
    stacked_data <- rbind(stacked_data,
                          data.frame(Comparison = comp_name,
                                     Cluster = new_cluster,
                                     sig_count = count_sig,
                                     stringsAsFactors = FALSE))
  }
}

# For each comparison, compute the percentage contribution of each cluster.
stacked_data <- stacked_data %>% 
  group_by(Comparison) %>%
  mutate(Percentage = sig_count / sum(sig_count) * 100)

# Reorder clusters so that the clusters with lower overall significant counts appear on top.
# Here, we calculate the total sig_count across all comparisons, sort in ascending order,
# and then reverse the order to have the lowest percentage plotted last (on top of the stack).
cluster_order <- stacked_data %>% 
  group_by(Cluster) %>% 
  summarise(total = sum(sig_count)) %>% 
  arrange(total) %>% 
  pull(Cluster)
stacked_data$Cluster <- factor(stacked_data$Cluster, levels = rev(cluster_order))

# Create the stacked bar plot.
stacked_plot <- ggplot(stacked_data, aes(x = Comparison, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Distribution of Significant Exons by Cluster",
       y = "Percentage of Significant Events (%)",
       x = "Comparison") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

# Display the plot.
print(stacked_plot)



# Initialize an empty list to store the plots.
plots <- list()

# Load necessary packages
library(ggplot2)
library(dplyr)
library(patchwork)
library(showtext)

font_add_google("Courier Prime", "courier")
showtext_auto()
showtext_opts(dpi=300)

# Your color/alpha palettes
color_palette <- c("FALSE" = "#9E9E9E", "TRUE"  = "#00BFFF")
alpha_palette <- c("FALSE" = 0.2,      "TRUE"  = 1)

ncol <- 2
bottom_row <- (length(clusters) - ncol + 1):length(clusters)

plots <- vector("list", length(clusters))

nrow <- ceiling(length(clusters) / ncol)  # Number of rows in the final grid

for(i in seq_along(clusters)) {
  cl         <- clusters[i]
  gene_list  <- clust_exons[[cl]]
  event_count<- sub(".*\\((\\d+) genes\\).*", "\\1", cl)
  new_title  <- paste0("Cluster ", i, " (", event_count, " Events)")
  
  p <- ggplot(
    filter(exons_longer, Genes %in% gene_list),
    aes(x = as.factor(time), y = value, group = Genes, color = significant)
  ) +
    geom_line(aes(alpha = significant), size = 0.8) +
    scale_color_manual(values = color_palette) +
    scale_alpha_manual(values = alpha_palette) +
    scale_x_discrete(
      breaks = c("t0", "t2", "t4", "t7"),
      labels = c("D0", "D2", "D4", "D7")
    ) +
    labs(
      title = new_title,
      y     = "Normalised PSI",
      x     = "Days of Differentiation"
    ) +
    theme(
      text              = element_text(family = "courier"),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black", size = 0.3),
      axis.ticks        = element_line(size = 1.5),
      axis.title.y      = element_text(size = 18, face = "bold"),
      axis.text.y       = element_text(size = 18),
      axis.text.x       = element_text(size = 18),
      title = element_text(size=20),
      legend.position   = "none"
    )
  
  # Remove x-axis title if not in the last row
  if (i <= (nrow - 1) * ncol) {
    p <- p + theme(axis.title.x = element_blank())
  } else {
    p <- p + theme(axis.title.x = element_text(size = 18, face = "bold"))
  }
  
  plots[[i]] <- p
}

# Combine with patchwork
composite_plot_exons <- wrap_plots(plots, ncol = ncol) +
  plot_annotation(
    theme = theme(
      text            = element_text(family = "courier"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

print(composite_plot_exons)


ggsave("composite_plot_exons.png",
       plot = composite_plot_exons,
       width = 14, height =12 , dpi = 300)


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

plots <- vector("list", length(clusters))

nrow <- ceiling(length(clusters) / ncol)  # Number of rows in the final grid

for(i in seq_along(clusters)) {
  cl         <- clusters[i]
  gene_list  <- clust_introns[[cl]]
  event_count<- sub(".*\\((\\d+) genes\\).*", "\\1", cl)
  new_title  <- paste0("Cluster ", i, " (", event_count, " Events)")
  
  p <- ggplot(
    filter(introns_longer, Genes %in% gene_list),
    aes(x = as.factor(time), y = value, group = Genes, color = significant)
  ) +
    geom_line(aes(alpha = significant), size = 0.8) +
    scale_color_manual(values = color_palette) +
    scale_alpha_manual(values = alpha_palette) +
    scale_x_discrete(
      breaks = c("t0", "t2", "t4", "t7"),
      labels = c("D0", "D2", "D4", "D7")
    ) +
    labs(
      title = new_title,
      y     = "Normalised PSI",
      x     = "Days of Differentiation"
    ) +
    theme(
      text              = element_text(family = "courier"),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black", size = 0.3),
      axis.ticks        = element_line(size = 1.5),
      axis.title.y      = element_text(size = 18, face = "bold"),
      axis.text.y       = element_text(size = 18),
      axis.text.x       = element_text(size = 18),
      title = element_text(size=20),
      legend.position   = "none"
    )
  
  # Remove x-axis title if not in the last row
  if (i <= (nrow - 1) * ncol) {
    p <- p + theme(axis.title.x = element_blank())
  } else {
    p <- p + theme(axis.title.x = element_text(size = 18, face = "bold"))
  }
  
  plots[[i]] <- p
}

# Combine with patchwork
composite_plot_introns <- wrap_plots(plots, ncol = ncol) +
  plot_annotation(
    theme = theme(
      text            = element_text(family = "courier"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

print(composite_plot_introns)


ggsave("composite_plot_introns.png",
       plot = composite_plot_introns,
       width = 14, height =6 , dpi = 300)



