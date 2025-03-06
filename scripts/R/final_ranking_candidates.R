#####
# This code coverts the final list fo candaidtes with the three flags to a ranked database based on:
# 1) DeltaPSI -- from 0 to 1 score
# 2) pvalue -- from 0 to 1 score using the log conversion
# 3) Number of hits -- from 0 to (number of hits - 1)
# 4) Whether the exon has inclusion (dPSI>0) -- 0 or 1 if T
# 5) Follows a linear depndency from 0 days to 2 days to 4 days -- from 0 to 1 ranked [dPSI0->2 * dPSI 2->4]
# 6) It is highly expressed in C2C12 cells -- 0 or 1 if T
#####
final_list<-read.csv("final_good_list.csv")[,-1]

library(jsonlite)
library(tidyverse)
library(biomaRt)

# Load the JSON file
file_path <- "BioGPS+Mouse+Cell+Type+and+Tissue+Gene+Expression+Profiles.txt"  # Replace with the actual file path
json_data <- fromJSON(file_path)



tao_pdiff_exons_2d<-na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_exons_2d.csv")[,-1])
tao_pdiff_exons_4dvs2d<-na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_exons_4dvs2d.csv")[,-1])
linear_df<-data.frame(event= tao_pdiff_exons_4dvs2d$EVENT, linearity= tao_pdiff_exons_2d$deltapsi*tao_pdiff_exons_4dvs2d$deltapsi)
tao_pdiff_introns_2d<-na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_introns_2d.csv")[,-1])
tao_pdiff_introns_4dvs2d<-na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_introns_4dvs2d.csv")[,-1])
linear_df<-rbind(linear_df, data.frame(event= tao_pdiff_introns_4dvs2d$EVENT, linearity= tao_pdiff_introns_2d$deltapsi*tao_pdiff_introns_4dvs2d$deltapsi))

diff_express<-read.csv("DESeq2_d4_vs_d0.csv")
diff_express<-na.omit(diff_express)

highExpGenes<-readRDS("highExpGenes.rds")
lowExpGenes<-readRDS("lowExpGenes.rds")

ranking_df <- final_list %>%
  mutate(ranked_deltapsi = abs(deltapsi) / max(abs(deltapsi))) %>%
  mutate(exon_included = ifelse(deltapsi > 0, 1, 0)) %>%
  # Get the corresponding linearity value from linear_df
  mutate(raw_linearity = linear_df$linearity[match(EVENT, linear_df$event)]) %>%
  # For negatives (or zero) assign 0; for positives, take a log-transform (log(x+1))
  mutate(linear_dependent = ifelse(raw_linearity > 0,
                                   log(raw_linearity + 1),
                                   0)) %>%
  # Rescale the positive (nonzero) log-transformed values to be between 0 and 1
  mutate(linear_dependent = ifelse(linear_dependent > 0,
                                   (linear_dependent - min(linear_dependent[linear_dependent > 0], na.rm = TRUE)) /
                                     (max(linear_dependent[linear_dependent > 0], na.rm = TRUE) - min(linear_dependent[linear_dependent > 0], na.rm = TRUE)),
                                   0)) %>%
  mutate(expression = diff_express$log2FoldChange[match(GENE, diff_express$Gene)]) %>%
  mutate(more_expressed = ifelse((abs(expression) < 0.5) & (GENE %in% highExpGenes), 1,
                                 ifelse(abs(expression) >= 0.5 & (GENE %in% lowExpGenes), -1, 0))) %>%
  group_by(GENE) %>%
  mutate(hits = n()) %>%
  ungroup() %>%
  select(GENE, EVENT, LENGTH, sequence, ranked_deltapsi, exon_included, linear_dependent, more_expressed, hits)


# Compute scaled variables and composite score
score_df <- ranking_df %>%
  mutate(
    more_expressed = (more_expressed + 1) / 2,
    hits_scaled = hits / 4,                        # Scale hits from [1,4] to [0.25,1]
    composite_score = ranked_deltapsi + exon_included + 
      linear_dependent + hits_scaled + more_expressed,
    candidate = paste(GENE, EVENT, sep = "_")
  ) %>%
  group_by(candidate) %>% 
  summarise(across(c(ranked_deltapsi, exon_included, linear_dependent, 
                    hits_scaled, composite_score, more_expressed), mean)) %>%
  arrange(composite_score)

# Reshape to long format for a stacked bar chart
long_df <- score_df %>%
  select(candidate, ranked_deltapsi, exon_included, linear_dependent, hits_scaled, more_expressed) %>%
  pivot_longer(cols = -candidate, names_to = "metric", values_to = "value")

# Rename factors for clarity and ordering in the legend
long_df$metric <- factor(long_df$metric,
                         levels = c("hits_scaled", "more_expressed","linear_dependent", "exon_included", "ranked_deltapsi"),
                         labels = c("Hits (scaled)", "Higher Expression (4 days vs 0 Days)", "Linear Dependent", "Exon Included", "Ranked DeltaPSI"))

# Reorder candidate factor based on composite score
long_df$candidate <- factor(long_df$candidate, levels = unique(score_df$candidate))

# Define custom colors for each metric
metric_colors <- c("Hits (scaled)" = "#66c2a5", 
                   #"High Absolute Expression (C2C12)" = "#fc8d62",
                   "Higher Expression (4 days vs 0 Days)" = "#f1c95f",
                   "Linear Dependent" = "#8da0cb", 
                   "Exon Included" = "#e78ac3", 
                   "Ranked DeltaPSI" = "#a6d854")

# Create the stacked bar plot
p <- ggplot(long_df, aes(x = candidate, y = value, fill = metric)) +
  geom_bar(stat = "identity") +
  # Add total composite score as text labels correctly
  geom_text(data = score_df, aes(x = candidate, y = composite_score + 0.05, 
                                 label = round(composite_score, 2)), inherit.aes = FALSE,
            hjust = -0.1, size = 3, color = "black") +
  geom_hline(yintercept = 6, linewidth = 0.8) +
  # Add label above the hline
  geom_text(aes(x = 1.3, y = 5.3, label = "Max Theoretical Value"), color = "black", 
            size = 4, fontface = "italic", hjust = 0) +
  coord_flip() +  # Horizontal bars for better readability
  scale_fill_manual(values = metric_colors) +
  labs(title = "Composite Scores of Candidates",
       subtitle = "Stacked contributions from Ranked DeltaPSI, Exon Inclusion, Linear Dependency, High Expression, and Hits",
       x = "Candidate (GENE_EVENT)",
       y = "Contribution to Composite Score",
       fill = "Metric") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.position = "top",
        plot.title = element_text(face = "bold")) +
  # Extend y-axis limit for clear text display
  ylim(0, max(score_df$composite_score) * 1.1)


# Display the plot
print(p)

