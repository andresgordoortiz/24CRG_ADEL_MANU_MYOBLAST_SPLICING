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

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Read and process exon data
tao_pdiff_exons_2d <- na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_exons_2d.csv")[, -1])
tao_pdiff_exons_4dvs2d <- na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_exons_4dvs2d.csv")[, -1])
linear_df <- data.frame(
  event = tao_pdiff_exons_4dvs2d$EVENT, 
  linearity = tao_pdiff_exons_2d$deltapsi * tao_pdiff_exons_4dvs2d$deltapsi
)

# Read and process intron data and combine with exon linearity values
tao_pdiff_introns_2d <- na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_introns_2d.csv")[, -1])
tao_pdiff_introns_4dvs2d <- na.omit(read.csv("notebooks/final/tables_10n/tao_pdiff_introns_4dvs2d.csv")[, -1])
linear_df <- rbind(
  linear_df,
  data.frame(
    event = tao_pdiff_introns_4dvs2d$EVENT, 
    linearity = tao_pdiff_introns_2d$deltapsi * tao_pdiff_introns_4dvs2d$deltapsi
  )
)

# Read differential expression and gene expression data
diff_express <- na.omit(read.csv("DESeq2_d4_vs_d0.csv"))
highExpGenes <- readRDS("highExpGenes.rds")
lowExpGenes <- readRDS("lowExpGenes.rds")

# Process the main ranking data (assumes final_list exists in the environment)
ranking_df <- final_list %>%
  mutate(
    ranked_deltapsi = abs(deltapsi) / max(abs(deltapsi)),
    exon_included = ifelse(deltapsi > 0, 1, 0),
    # Look up the corresponding linearity value from linear_df
    raw_linearity = linear_df$linearity[match(EVENT, linear_df$event)],
    # Apply log-transform for positive values; otherwise assign 0
    linear_dependent = ifelse(raw_linearity > 0, log(raw_linearity + 1), 0)
  ) %>%
  mutate(
    # Rescale positive (nonzero) log-transformed values to [0, 1]
    linear_dependent = ifelse(
      linear_dependent > 0,
      (linear_dependent - min(linear_dependent[linear_dependent > 0], na.rm = TRUE)) /
        (max(linear_dependent[linear_dependent > 0], na.rm = TRUE) -
           min(linear_dependent[linear_dependent > 0], na.rm = TRUE)),
      0
    ),
    # Match gene expression values
    expression = diff_express$log2FoldChange[match(GENE, diff_express$Gene)],
    more_expressed = ifelse(abs(expression) < 0.5, 1, 0),
    high_expressed = ifelse(GENE %in% highExpGenes, 1, ifelse(GENE %in% lowExpGenes, -1, 0)),
    # Compute length rank differently for exons vs. introns
    length_rank = case_when(
      grepl("EX", EVENT) ~ (LENGTH - min(LENGTH[grepl("EX", EVENT)])) /
        (max(LENGTH[grepl("EX", EVENT)]) - min(LENGTH[grepl("EX", EVENT)])),
      grepl("INT", EVENT) ~ 1 - ((LENGTH - min(LENGTH[grepl("INT", EVENT)])) /
                                   (max(LENGTH[grepl("INT", EVENT)]) - min(LENGTH[grepl("INT", EVENT)]))),
      TRUE ~ (LENGTH - min(LENGTH)) / (max(LENGTH) - min(LENGTH))
    )
  ) %>%
  group_by(GENE) %>%
  mutate(hits = n()) %>%
  ungroup() %>%
  select(GENE, EVENT, LENGTH, sequence, ranked_deltapsi, exon_included,
         linear_dependent, more_expressed, hits, length_rank, high_expressed)

# Compute scaled variables and composite score
score_df <- ranking_df %>%
  mutate(
    hits_scaled = hits / max(hits),
    high_expressed = (high_expressed + 1) / 2, 
    composite_score = ranked_deltapsi + exon_included + 
      linear_dependent + hits_scaled + more_expressed + length_rank + high_expressed,
    candidate = paste(GENE, EVENT, sep = "_")
  ) %>%
  group_by(candidate) %>% 
  summarise(across(
    c(ranked_deltapsi, exon_included, linear_dependent, hits_scaled, 
      composite_score, more_expressed, length_rank, high_expressed), mean
  )) %>%
  arrange(composite_score)

# Reshape to long format for a stacked bar chart
long_df <- score_df %>%
  select(candidate, ranked_deltapsi, exon_included, linear_dependent,
         hits_scaled, more_expressed, length_rank, high_expressed) %>%
  pivot_longer(cols = -candidate, names_to = "metric", values_to = "value")

# Rename factors for clarity and ordering in the legend
long_df$metric <- factor(
  long_df$metric,
  levels = c("hits_scaled", "length_rank", "more_expressed", "high_expressed",
             "linear_dependent", "exon_included", "ranked_deltapsi"),
  labels = c(
    "Number of Hits (scaled)",
    "Length Ranked 0->1 (Exons and Introns Separately)",
    "Difference in Expression 0:Yes, 1:No",
    "High Abundance 0:Q25, 0.5:Q25-Q75, 1:Q75",
    "Linear Dependent Ranked 0->1 (log Scale)",
    "ΔPSI > 0 (Inclusion) 0:No, 1:Yes",
    "ΔPSI Ranked 0->1"
  )
)

# Reorder candidate factor based on composite score
long_df$candidate <- factor(long_df$candidate, levels = unique(score_df$candidate))

# Define custom colors for each metric (added a color for "High Abundance")
metric_colors <- c(
  "Number of Hits (scaled)" = "#66c2a5", 
  "Length Ranked 0->1 (Exons and Introns Separately)" = "#fc8d62",
  "Difference in Expression 0:Yes, 1:No" = "#f1c95f",
  "High Abundance 0:Q25, 0.5:Q25-Q75, 1:Q75" = "#D95F02",
  "Linear Dependent Ranked 0->1 (log Scale)" = "#8da0cb", 
  "ΔPSI > 0 (Inclusion) 0:No, 1:Yes" = "#e78ac3", 
  "ΔPSI Ranked 0->1" = "#a6d854"
)

# Create the stacked bar plot with corrected vline for the theoretical max
p <- ggplot(long_df, aes(x = candidate, y = value, fill = metric)) +
  geom_bar(stat = "identity") +
  # Add composite score labels
  geom_text(
    data = score_df,
    aes(x = candidate, y = composite_score + 0.05, label = round(composite_score, 2)),
    inherit.aes = FALSE, hjust = -0.1, size = 3, color = "black"
  ) +
  # Add a vertical line at x = 7 (the theoretical max composite score)
  geom_hline(yintercept = 7, linewidth = 0.8, color = "black") +
  # Annotate the vertical line
  annotate(
    "text",
    x = length(levels(score_df$candidate)) + 0.5,
    y = 6.5,  # positions text above the candidate list
    label = "Max Theoretical Value",
    vjust = -0.5,
    size = 5,
    fontface = "italic",
    color = "black"
  ) +
  coord_flip() +  # Flip coordinates so candidate is on the vertical axis
  scale_fill_manual(values = metric_colors) +
  labs(
    title = "Composite Scores of Candidates",
    subtitle = "Stacked contributions from Ranked DeltaPSI, Exon Inclusion, Linear Dependency, High Expression, and Hits",
    x = "Candidate (GENE_EVENT)",
    y = "Contribution to Composite Score",
    fill = "Metric"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) 

# Display the plot
print(p)

