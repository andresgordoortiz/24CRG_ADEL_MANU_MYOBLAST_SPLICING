tao_data <- getDataset(pathTables = paste0(getwd(),"/inclusion_tables/Tao_INCLUSION_LEVELS_FULL-mm10.tab"), tool = "vast-tools")
tao_events <- filterEvents(getEvents(tao_data, tool = "vast-tools"), N=10) # Extract alternative splicing events
tao_exons <- filterEvents(tao_events, types = c("C1", "C2", "C3", "S", "MIC"), N = 10)
tao_introns <- filterEvents(tao_events, types = c("IR"), N = 10)

splice_factors<-read_xlsx("SFs_list_ensembl_ID.xlsx")

correlation_matrix<-tao_events$PSI[,-c(3:6)]



correlation_matrix<-correlation_matrix[, c(colnames(correlation_matrix)[1:2], metadata_tao$run_accession)]
subset_genes <- correlation_matrix[tolower(correlation_matrix$GENE) %in% tolower(splice_factors$symbol), ]
rest_genes <- correlation_matrix[!tolower(correlation_matrix$GENE) %in% tolower(splice_factors$symbol), ]
subset_matrix <- as.matrix(subset_genes[, 3:14])
rest_matrix <- as.matrix(rest_genes[, 3:14])


results <- data.frame(
  Gene1 = character(),
  Gene2 = character(),
  Correlation = numeric(),
  P_value = numeric(),
  Adjusted_P_value = numeric(),
  stringsAsFactors = FALSE
)

# Compute correlations and p-values
all_p_values <- c()
for (i in 1:nrow(subset_matrix)) {
  for (j in 1:nrow(rest_matrix)) {
    # Perform correlation test
    test <- cor.test(as.numeric(subset_matrix[i, ]), as.numeric(rest_matrix[j, ]), method = "pearson", use = "complete.obs")
    
    # Store the raw p-values for later adjustment
    all_p_values <- c(all_p_values, test$p.value)
    
    # Temporarily store the results (to be adjusted later)
    results <- rbind(results, data.frame(
      Gene1 = subset_genes$GENE[i],
      Gene2 = rest_genes$GENE[j],
      Correlation = test$estimate,
      P_value = test$p.value,
      Adjusted_P_value = NA  # Placeholder for the adjusted p-value
    ))
  }
}

# Apply multiple testing correction (Benjamini-Hochberg)
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

# Update the results data frame with adjusted p-values
results$Adjusted_P_value <- adjusted_p_values

# Filter significant results (adjusted p-value < 0.05)
significant_results <- na.omit(results[results$Adjusted_P_value < 0.01 & abs(results$Correlation) >=0.9, ])

write.csv(significant_results, "significant_results.csv")
