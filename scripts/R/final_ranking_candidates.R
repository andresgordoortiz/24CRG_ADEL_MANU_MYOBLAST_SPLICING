#####
# This code coverts the final list fo candaidtes with the three flags to a ranked database based on:
# 1) DeltaPSI -- from 0 to 1 score
# 2) pvalue -- from 0 to 1 score using the log conversion
# 3) Number of hits -- from 0 to (number of hits - 1)
# 4) Whether the exon has inclusion (dPSI>0) -- 0 or 1 if T
# 5) Follows a linear depndency from 0 days to 2 days to 4 days -- from 0 to 1 ranked [dPSI0->2 * dPSI 2->4]
# 6) It is highly expressed in C2C12 cells -- 0 or 1 if T
#####

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


##### 
# All Layers
#####
final_all_layers<-read.csv("final_all_layers.csv")[,-1]
all_layers<-final_all_layers
# Process the main ranking data (assumes final_list exists in the environment)
ranking_df <- final_all_layers %>%
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
p_all <- ggplot(long_df, aes(x = candidate, y = value, fill = metric)) +
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
    title = "Composite Scores of Candidates in All Layers",
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
print(p_all)

#####
# Explanation: 
# *Number of Hits*: Number of spliced events in the candidate gene. The higher the number, the more differentially spliced a given gene is.
# Length Ranked: RNA FISH requires long enough target sequences for probe design. Thus, longer exons will be given priority.
# Difference in Expression: Whether or not the gene is differentially (log2FC<0.5) expressed comparing day 4 to day 0 of myogenesis.
# High Abundance: Whether the gene is in the top 25% of expression (0), between 25% and 75% (0.5), or in the top 25% (1).
# Linear Dependent: The linearity of the splicing pattern across the time points (day 2 to day 0, and day 4 to day 2). The more linear, the higher the score.
# DeltaPSI >0: Whether the gene is included (1) or excluded (0) in the mature mRNA. Nuclear agitation seems to promote exon inclusion.
# deltaPSI ranked: The absolute value of the deltaPSI, normalized to the maximum value. The higher the value, the more differentially spliced the gene is.
#####

##### 
# Myogenesis & Forces
#####
final_myo_forces_layers<-read.csv("final_myo_forces_layers.csv")[,-1]
myo_forces_layers<-final_myo_forces_layers
# Process the main ranking data (assumes final_list exists in the environment)
ranking_df <- final_myo_forces_layers %>%
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
    more_expressed = ifelse(is.na(expression), 0, 
                            ifelse(abs(expression) < 0.5, 1, 0)),
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

# Filter only the top 20 candidates
top_n <- 30
filtered_long_df <- long_df %>% filter(candidate %in% tail(unique(candidate), top_n))
filtered_score_df <- score_df %>% filter(candidate %in% tail(unique(candidate), top_n))

# Count omitted candidates
omitted_count <- length(unique(long_df$candidate)) - top_n

# Create the stacked bar plot with corrected vline for the theoretical max
p_myo_forces <- ggplot(filtered_long_df, aes(x = candidate, y = value, fill = metric)) +
  geom_bar(stat = "identity") +
  # Add composite score labels
  geom_text(
    data = filtered_score_df,
    aes(x = candidate, y = composite_score + 0.05, label = round(composite_score, 2)),
    inherit.aes = FALSE, hjust = -0.1, size = 3, color = "black"
  ) +
  # Add a vertical line at x = 7 (the theoretical max composite score)
  geom_hline(yintercept = 7, linewidth = 0.8, color = "black") +
  # Annotate the vertical line
  annotate(
    "text",
    x = 1,  # Adjusted for filtered candidates
    y = 6.5,  # positions text above the candidate list
    label = "Max Theoretical Value",
    vjust = -0.5,
    size = 5,
    fontface = "italic",
    color = "black"
  ) +
  # Annotate the number of omitted candidates
  annotate(
    "text",
    x = top_n,
    y = max(filtered_score_df$composite_score, na.rm = TRUE) + 0.8,
    label = paste0("(", omitted_count, " candidates omitted)"),
    size = 4,
    color = "gray50",
    fontface = "italic"
  ) +
  coord_flip() +  # Flip coordinates so candidate is on the vertical axis
  scale_fill_manual(values = metric_colors) +
  labs(
    title = "Composite Scores of Candidates in Myogenesis & Force-affected",
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

# Print the plot
p_myo_forces


##### 
# Myogenesis & Speckle
#####
final_myo_associated_layers<-read.csv("final_myo_associated_layers.csv")[,-1]
myo_associated_layers<-final_myo_associated_layers
# Process the main ranking data (assumes final_list exists in the environment)
ranking_df <- final_myo_associated_layers %>%
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
    more_expressed = ifelse(is.na(expression), 0, 
                            ifelse(abs(expression) < 0.5, 1, 0)),
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

# Filter only the top 20 candidates
top_n <- 30
filtered_long_df <- long_df %>% filter(candidate %in% tail(unique(candidate), top_n))
filtered_score_df <- score_df %>% filter(candidate %in% tail(unique(candidate), top_n))

# Count omitted candidates
omitted_count <- length(unique(long_df$candidate)) - top_n

# Create the stacked bar plot with corrected vline for the theoretical max
p_myo_speckle <- ggplot(filtered_long_df, aes(x = candidate, y = value, fill = metric)) +
  geom_bar(stat = "identity") +
  # Add composite score labels
  geom_text(
    data = filtered_score_df,
    aes(x = candidate, y = composite_score + 0.05, label = round(composite_score, 2)),
    inherit.aes = FALSE, hjust = -0.1, size = 3, color = "black"
  ) +
  # Add a vertical line at x = 7 (the theoretical max composite score)
  geom_hline(yintercept = 7, linewidth = 0.8, color = "black") +
  # Annotate the vertical line
  annotate(
    "text",
    x = 1,  # Adjusted for filtered candidates
    y = 6.5,  # positions text above the candidate list
    label = "Max Theoretical Value",
    vjust = -0.5,
    size = 5,
    fontface = "italic",
    color = "black"
  ) +
  # Annotate the number of omitted candidates
  annotate(
    "text",
    x = top_n,
    y = max(filtered_score_df$composite_score, na.rm = TRUE) + 0.8,
    label = paste0("(", omitted_count, " candidates omitted)"),
    size = 4,
    color = "gray50",
    fontface = "italic"
  ) +
  coord_flip() +  # Flip coordinates so candidate is on the vertical axis
  scale_fill_manual(values = metric_colors) +
  labs(
    title = "Composite Scores of Candidates in Myogenesis & Speckle-Associated",
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

# Print the plot
p_myo_speckle

##### 
# Myogenesis Only
#####
final_only_myo_layers<-read.csv("final_only_myo_layers.csv")[,-1]
only_myo_layers<-final_only_myo_layers
# Process the main ranking data (assumes final_list exists in the environment)
ranking_df <- final_only_myo_layers %>%
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
    more_expressed = ifelse(is.na(expression), 0, 
                            ifelse(abs(expression) < 0.5, 1, 0)),
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

# Filter only the top 20 candidates
top_n <- 30
filtered_long_df <- long_df %>% filter(candidate %in% tail(unique(candidate), top_n))
filtered_score_df <- score_df %>% filter(candidate %in% tail(unique(candidate), top_n))

# Count omitted candidates
omitted_count <- length(unique(long_df$candidate)) - top_n

# Create the stacked bar plot with corrected vline for the theoretical max
p_myo_only <- ggplot(filtered_long_df, aes(x = candidate, y = value, fill = metric)) +
  geom_bar(stat = "identity") +
  # Add composite score labels
  geom_text(
    data = filtered_score_df,
    aes(x = candidate, y = composite_score + 0.05, label = round(composite_score, 2)),
    inherit.aes = FALSE, hjust = -0.1, size = 3, color = "black"
  ) +
  # Add a vertical line at x = 7 (the theoretical max composite score)
  geom_hline(yintercept = 7, linewidth = 0.8, color = "black") +
  # Annotate the vertical line
  annotate(
    "text",
    x = 1,  # Adjusted for filtered candidates
    y = 6.5,  # positions text above the candidate list
    label = "Max Theoretical Value",
    vjust = -0.5,
    size = 5,
    fontface = "italic",
    color = "black"
  ) +
  # Annotate the number of omitted candidates
  annotate(
    "text",
    x = top_n,
    y = max(filtered_score_df$composite_score, na.rm = TRUE) + 0.8,
    label = paste0("(", omitted_count, " candidates omitted)"),
    size = 4,
    color = "gray50",
    fontface = "italic"
  ) +
  coord_flip() +  # Flip coordinates so candidate is on the vertical axis
  scale_fill_manual(values = metric_colors) +
  labs(
    title = "Composite Scores of Candidates in Myogenesis Only",
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

# Print the plot
p_myo_only


#####
# Exploratory Table
#####

# ---------------------------
# 1. Load required libraries
# ---------------------------
library(dplyr)
library(knitr)
library(ggplot2)
library(gridExtra)

# ---------------------------
# 2. Define helper functions
# ---------------------------

# Function to calculate GC content (in percentage) from a DNA sequence
calc_gc <- function(seq) {
  seq <- toupper(seq)
  nucleotides <- unlist(strsplit(seq, split = ""))
  gc_count <- sum(nucleotides %in% c("G", "C"))
  return(100 * gc_count / nchar(seq))
}

# Function to summarize a data frame for rows that match an event pattern ("EX" or "INT")
summarize_by_event <- function(df, event_pattern) {
  df_subset <- df %>% filter(grepl(event_pattern, EVENT))
  
  if(nrow(df_subset) > 0) {
    gc_values <- sapply(df_subset$sequence, calc_gc)
    summary_stats <- df_subset %>% summarise(
      LENGTH_mean   = mean(LENGTH, na.rm = TRUE),
      LENGTH_median = median(LENGTH, na.rm = TRUE),
      LENGTH_sd     = sd(LENGTH, na.rm = TRUE),
      GC_mean       = mean(gc_values, na.rm = TRUE),
      GC_median     = median(gc_values, na.rm = TRUE),
      GC_sd         = sd(gc_values, na.rm = TRUE)
    )
  } else {
    summary_stats <- data.frame(
      LENGTH_mean   = NA,
      LENGTH_median = NA,
      LENGTH_sd     = NA,
      GC_mean       = NA,
      GC_median     = NA,
      GC_sd         = NA
    )
  }
  
  return(summary_stats)
}

# ---------------------------
# 3. Compute Summary Statistics
# ---------------------------

# (Assuming your four data frames already exist: all_layers, myo_forces_layers, myo_associated_layers, only_myo_layers)

# Create summary for EX events for each data frame
summary_EX <- bind_rows(
  summarize_by_event(all_layers, "EX") %>% mutate(DataFrame = "All_Layers", Event_Type = "EX"),
  summarize_by_event(myo_forces_layers, "EX") %>% mutate(DataFrame = "Myo_Forces_Layers", Event_Type = "EX"),
  summarize_by_event(myo_associated_layers, "EX") %>% mutate(DataFrame = "Myo_Associated_Layers", Event_Type = "EX"),
  summarize_by_event(only_myo_layers, "EX") %>% mutate(DataFrame = "Only_Myo_Layers", Event_Type = "EX")
)

# Create summary for INT events for each data frame
summary_INT <- bind_rows(
  summarize_by_event(all_layers, "INT") %>% mutate(DataFrame = "All_Layers", Event_Type = "INT"),
  summarize_by_event(myo_forces_layers, "INT") %>% mutate(DataFrame = "Myo_Forces_Layers", Event_Type = "INT"),
  summarize_by_event(myo_associated_layers, "INT") %>% mutate(DataFrame = "Myo_Associated_Layers", Event_Type = "INT"),
  summarize_by_event(only_myo_layers, "INT") %>% mutate(DataFrame = "Only_Myo_Layers", Event_Type = "INT")
)

# Display the summary tables using knitr::kable
kable(summary_EX, caption = "EX Summary Statistics")
kable(summary_INT, caption = "INT Summary Statistics")

# ---------------------------
# 4. Prepare Data for ANOVA & Plotting
# ---------------------------

# Create a named list of your four data frames
df_list <- list(
  All_Combined = all_layers,
  Myogenesis_Forces = myo_forces_layers,
  Myogenesis_Speckles = myo_associated_layers,
  Only_Myogenesis = only_myo_layers
)

# Function to combine data for a given event type for plotting
combine_event_data <- function(df_list, event_pattern) {
  combined_data <- do.call(rbind, lapply(names(df_list), function(name){
    df <- df_list[[name]]
    
    # Ensure the required columns exist
    if (!all(c("LENGTH", "EVENT", "sequence") %in% colnames(df))) {
      return(NULL)
    }
    
    # Filter for the event type (EX or INT)
    df_sub <- df %>% filter(grepl(event_pattern, EVENT))
    
    if(nrow(df_sub) > 0) {
      df_sub <- df_sub %>%
        mutate(group = name,  # Assign the data frame name as a group label
               GC = sapply(sequence, calc_gc)) %>% 
        select(LENGTH, GC, group) # Keep only relevant columns
      return(df_sub)
    } else {
      return(NULL)
    }
  }))
  return(combined_data)
}


# Combine data for EX and INT events
combined_ex <- combine_event_data(df_list, "EX")
combined_int <- combine_event_data(df_list, "INT")

# ---------------------------
# 5. Run ANOVA Tests
# ---------------------------

# ANOVA for EX events
# Note that variances are not equal, thus the ANOVA Welch
anova_ex_length <- oneway.test(LENGTH ~ group, data = combined_ex,var.equal = FALSE)[[3]]
anova_ex_gc     <- oneway.test(GC ~ group, data = combined_ex,var.equal = FALSE)[[3]]

# ANOVA for INT events
anova_int_length <- oneway.test(LENGTH ~ group, data = combined_int, var.equal = FALSE)[[3]]
anova_int_gc     <- oneway.test(GC ~ group, data = combined_int, var.equal = FALSE)[[3]]

# ---------------------------
# 6. Create Boxplots for ANOVA Visualization
# ---------------------------

library(ggplot2)
library(ggpubr)
library(extrafont)
# Ensure fonts are imported and loaded
# font_import(pattern = "Arial", prompt = FALSE) # Run once if needed
loadfonts(device = "win")

# Define a refined theme for publication-quality plots
theme_pub <- theme_minimal(base_family = "sans") + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "#333333"),
    axis.title = element_text(size = 14, face = "bold", color = "#333333"),
    axis.text = element_text(size = 12, color = "#333333"),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Enhanced plotting function with GC-specific y-axis breaks and removal of non-significant comparison bars
create_boxplot <- function(data, x, y, title, p_value, y_transform = "none", y_limits = NULL) {
  # Define the specific pairwise comparisons
  comparisons_list <- list(
    c("All_Combined", "Myogenesis_Speckles"), 
    c("All_Combined", "Myogenesis_Forces"),
    c("All_Combined", "Only_Myogenesis"),
    c("Myogenesis_Speckles", "Myogenesis_Forces"),
    c("Myogenesis_Speckles", "Only_Myogenesis"),
    c("Myogenesis_Forces", "Only_Myogenesis")
  )
  
  # Compute pairwise comparisons using t.test
  comp_res <- compare_means(
    formula = as.formula(paste(y, "~", x)), 
    data = data,
    method = "t.test", 
    paired = FALSE,
    comparisons = comparisons_list
  )
  # Filter to include only significant comparisons (p < 0.05)
  signif_comps <- comp_res[comp_res$p < 0.05, c("group1", "group2")]
  if (nrow(signif_comps) > 0) {
    comparisons_to_plot <- apply(signif_comps, 1, function(x) c(x[1], x[2]))
    comparisons_to_plot <- split(t(comparisons_to_plot), seq(nrow(signif_comps)))
  } else {
    comparisons_to_plot <- NULL
  }
  
  # Build the plot
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, color = "black") + 
    geom_jitter(width = 0.2, alpha = 0.4, size = 2.5, color = "#444444") +
    # Add mean points for additional clarity
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
    scale_x_discrete(limits = c("All_Combined", "Myogenesis_Speckles", "Myogenesis_Forces", "Only_Myogenesis")) +
    labs(
      title = paste0(title, "\nANOVA Welch p-value = ", round(p_value, 4)),
      x = "Group",
      y = y
    ) +
    theme_pub
  
  # If the y variable is GC, set specific breaks to show percentages
  if (y == "GC") {
    p <- p + scale_y_continuous(breaks = c(0, 25, 50, 75, 100))
  }
  
  # Add pairwise comparisons only if there are significant ones
  if (!is.null(comparisons_to_plot)) {
    p <- p + stat_compare_means(
      comparisons = comparisons_to_plot, 
      method = "t.test", 
      label = "p.signif", 
      paired = FALSE, 
      var.equal = FALSE
    )
  }
  
  # Apply log10 transformation if specified (ensure y > 0)
  if (y_transform == "log10") {
    p <- p + scale_y_log10()
  }
  
  # Zoom in using specified y-axis limits if provided
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}

# Example usage for the LENGTH plots:
# Option 1: Use a log transformation (if appropriate)
p_length_ex <- create_boxplot(combined_ex, "group", "LENGTH", 
                              "LENGTH by Group (EX events)", anova_ex_length, 
                              y_transform = "log10")
p_length_int <- create_boxplot(combined_int, "group", "LENGTH", 
                               "LENGTH by Group (INT events)", anova_int_length, 
                               y_transform = "log10")

# Option 2: Zoom in on the bulk of the data using y-axis limits 
upper_limit_ex <- quantile(combined_ex$LENGTH, 0.99)
upper_limit_int <- quantile(combined_int$LENGTH, 0.99)

p_length_ex_zoom <- create_boxplot(combined_ex, "group", "LENGTH", 
                                   "LENGTH by Group (EX events)", anova_ex_length, 
                                   y_limits = c(min(combined_ex$LENGTH), upper_limit_ex))
p_length_int_zoom <- create_boxplot(combined_int, "group", "LENGTH", 
                                    "LENGTH by Group (INT events)", anova_int_length, 
                                    y_limits = c(min(combined_int$LENGTH), upper_limit_int))

# Generate the GC plots (the y-axis will now show only 0, 25, 50, 75, 100)
p_gc_ex  <- create_boxplot(combined_ex, "group", "GC", "GC Content by Group (EX events)", anova_ex_gc)
p_gc_int <- create_boxplot(combined_int, "group", "GC", "GC Content by Group (INT events)", anova_int_gc)

# Arrange the plots in a grid (choose the appropriate LENGTH plot version)
final_plot <- ggarrange(
  p_length_ex, p_gc_ex, p_length_int, p_gc_int, 
  ncol = 2, nrow = 2, 
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 16, face = "bold")
)

# Display the final arranged plot
final_plot




#####
# Inclusion patterns
######


# Function to compute the difference for a single dataframe
compute_difference <- function(df) {
  avg_4d <- rowMeans(df[, 9:11])
  avg_0d <- rowMeans(df[, 3:5])
  avg_4d - avg_0d
}



all_layers_matrix <- bind_rows(final_myo_forces_layers, final_all_layers) %>%
  distinct(.data$EVENT, .keep_all = TRUE) %>%
  select(1, 2, 7:18)


all_layers_matrix<-all_layers_matrix[, c("GENE","EVENT", metadata_tao$run_accession)]
# Calculate the difference and remove any additional NAs
all_layers_matrix$difference <- compute_difference(all_layers_matrix)
all_layers_matrix <- na.omit(all_layers_matrix)

# Create a group variable based on whether the difference is below or above zero
all_layers_matrix$group <- ifelse(all_layers_matrix$difference < 0, "Skipped", "Included")

# Perform a one-sample Wilcoxon signed-rank test (null hypothesis: median = 0)
symmetry_test <- wilcox.test(all_layers_matrix$difference, mu = 0)

# Compute the percentages (areas under the density curve)
pct_below <- round(mean(all_layers_matrix$difference < 0) * 100, 1)
pct_above <- round(mean(all_layers_matrix$difference > 0) * 100, 1)

# ---- Plot -------------------------------------------------------------------

plot_psi_distribution <- ggplot(
  all_layers_matrix[abs(all_layers_matrix$difference) > 0, ], 
  aes(x = difference, fill = group)
) +
  geom_histogram(aes(y = ..density..), bins = 200, 
                 position = "identity", alpha = 0.9) +
  labs(
    title = "Delta PSI Distribution of Significant Exons Myo and Forces",
    subtitle = paste("One-sample Wilcoxon test p-value for Symmetry:", 
                     format(symmetry_test$p.value, digits = 3)),
    x = "ΔPSI (4d - 0d)",
    y = "Density",
    fill = "Splicing Direction",
    caption = paste0("Created by AG on ", Sys.Date())
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "right",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(size = 15),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_x_continuous(
    breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100)
  ) +
  # Annotate the percentages with polished label boxes
  annotate("label", x = -65, y = 0.05, hjust = 0, vjust = 0,
           label = sprintf("dPSI < 0: %.1f%%", pct_below),
           fill = alpha("white", 0.8), color = "blue", size = 7, fontface = "bold",
           label.size = 0.5) +
  annotate("label", x = 65, y = 0.05, hjust = 1, vjust = 0,
           label = sprintf("dPSI > 0: %.1f%%", pct_above),
           fill = alpha("white", 0.8), color = "red", size = 7, fontface = "bold",
           label.size = 0.5) +
  scale_fill_brewer(palette = "Set1")

# Save and display the plot
ggsave("d4vs0_psi_distribution_exons.svg", plot_psi_distribution, width = 12, height = 8, dpi = 300)

print(plot_psi_distribution)


#####
# Gene centric Analysis
#####

# Set seed for reproducibility
set.seed(44)

# -------------------------------
# 0. Read protein impact data and combine your events
# -------------------------------
protein_impact <- read.table("notebooks/final/PROT_IMPACT-mm10-v2.3.tab", sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE, 
                             fill = TRUE, quote = "")

# Combine exon and intron event data (with a column "source" indicating the origin)
combined_df <- bind_rows(
  all  = final_all_layers,
  myo_speckle = final_myo_associated_layers,
  myo_force = final_myo_forces_layers,
  only_myo= final_only_myo_layers,
  .id = "source"
)

# Add protein impact info by matching on event IDs
combined_df$protein_impact <- protein_impact$ONTO[match(combined_df$EVENT, protein_impact$EventID)]

# Keep one row per EVENT (if duplicates exist)
combined_df <- combined_df %>% distinct(EVENT, .keep_all = TRUE) %>%
  filter(!GENE == "")

# -------------------------------
# 1. Create per-gene summary statistics
# -------------------------------
combined_df_summary <- combined_df %>%
  distinct(GENE, EVENT, .keep_all = TRUE) %>%  # remove duplicate gene/event rows if any
  group_by(GENE) %>%
  summarise(
    nEvents   = n(),                                        # number of events per gene
    meanDelta = mean(abs(deltapsi), na.rm = TRUE),          # mean absolute |dPSI|
    sdDelta   = sd(abs(deltapsi), na.rm = TRUE)             # standard deviation (for reference)
  ) %>%
  arrange(desc(nEvents), desc(meanDelta)) %>%              # sort: genes with more events first
  mutate(GENE = factor(GENE, levels = unique(GENE)))        # set ordering for later plotting

# -------------------------------
# 2. Add ORF information to gene-level summary
#    (We consider a gene “ORF disrupted” if any of its events has impact labeled as "ORF")
# -------------------------------
# First, mark each event as ORF or Non-ORF
combined_df <- combined_df %>%
  mutate(
    absDelta = abs(deltapsi),
    impact   = ifelse(!is.na(protein_impact) & grepl("ORF", protein_impact, ignore.case = TRUE),
                      "ORF", "Non-ORF")
  )

# First, compute per-gene ORF counts and the proportion of ORF events
gene_ORF_stats <- combined_df %>%
  distinct(GENE, EVENT, .keep_all = TRUE) %>%  # ensure one row per gene-event
  group_by(GENE) %>%
  summarise(
    total_events = n(),
    ORF_events   = sum(impact == "ORF", na.rm = TRUE),
    prop_ORF     = ORF_events / total_events
  )

# Merge these stats with your existing gene-level summary
combined_df_summary <- combined_df_summary %>%
  left_join(gene_ORF_stats, by = "GENE")

# Create the summary table, now using the proportion of ORF events
table_summary <- combined_df_summary %>%
  group_by(nEvents) %>%
  summarise(
    n_genes         = n(),
    avg_prop_ORF    = paste0(round(mean(prop_ORF, na.rm = TRUE), 3)*100, "%")
  ) %>%
  rename(
    "Number of Significant Events" = nEvents,
    "Number of Genes"              = n_genes,
    "Average Proportion of ORF"    = avg_prop_ORF
  )

library(DT)
datatable(table_summary, 
          rownames = FALSE,
          caption = 'Summary of Genes by Number of Significant Events')


# -------------------------------
# Filter for genes with >2 events and compute SD
# -------------------------------
df_summary_gt2 <- combined_df_summary %>%
  filter(nEvents > 2) %>%
  arrange(desc(nEvents), desc(meanDelta)) %>%
  mutate(
    GENE = factor(GENE, levels = unique(GENE))
  )

# -------------------------------
# Build the publication-ready plot:
# - Bars show the mean |dPSI| per gene.
# - Upper error bars show mean + SD.
# -------------------------------

p2 <- ggplot(df_summary_gt2, aes(x = GENE, y = meanDelta, fill = factor(nEvents))) +
  # Bar plot for mean |dPSI|
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  # Upper error bar: from meanDelta to meanDelta + SD
  geom_errorbar(aes(ymin = meanDelta, ymax = meanDelta + sdDelta),
                width = 0.2, color = "black", alpha = 0.6) +
  # Define fill colors (using a color palette) and label for the legend
  scale_fill_brewer(palette = "Set2", name = "Number of Events") +
  # Improved axis and plot annotations
  labs(
    x = "Gene ID",
    y = "Mean |dPSI|",
    title = "Genes with >2 Differential Splicing Events in Myoblasts",
    subtitle = "Bar height = mean |dPSI|; Upper error bar = SD"
  ) +
  # Use a clean, publication-ready theme with custom font settings
  theme_classic(base_size = 16, base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 20) # Ensure enough space
  )

# Display the plot
print(p2)
