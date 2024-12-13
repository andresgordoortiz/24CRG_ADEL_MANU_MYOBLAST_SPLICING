data<-rbind(zhang_pdiff_introns)
data<-data[abs(data$deltapsi)>=0.15,]
library(ggplot2)
library(dplyr)


# List of dataframes with metadata
dfs <- list(
  list(df = dominic_pdiff_exons, group = "Exons", study = "Dominic"),
  list(df = dominic_pdiff_introns, group = "Introns", study = "Dominic"),
  list(df = christopher_pdiff_exons, group = "Exons", study = "Christopher"),
  list(df = christopher_pdiff_introns, group = "Introns", study = "Christopher"),
  list(df = lyu_pdiff_exons, group = "Exons", study = "Lyu"),
  list(df = lyu_pdiff_introns, group = "Introns", study = "Lyu"),
  list(df = tao_pdiff_exons, group = "Exons", study = "Tao"),
  list(df = tao_pdiff_introns, group = "Introns", study = "Tao"),
  list(df = zhang_pdiff_exons, group = "Exons", study = "Zhang"),
  list(df = zhang_pdiff_introns, group = "Introns", study = "Zhang")
)

# Add metadata columns to each dataframe
dfs_with_metadata <- lapply(dfs, function(x) {
  x$df %>%
    mutate(Group = x$group, Study = x$study)
})

# Merge all dataframes by EVENT, filling missing columns with NA
final_df <- Reduce(function(x, y) {
  full_join(x, y, by = "EVENT")
}, dfs_with_metadata)

# View the final dataframe
print(final_df)



ggplot(data, aes(x=deltapsi, fill=ifelse(deltapsi < 0, "Skipped", "Included"))) +
  geom_histogram(alpha=0.7, bins=100) +
  theme_minimal() +
  labs(title="PSI distribution of Exons", x="ΔPSI (Myotubes-Myoblasts)", y="Density", fill="|ΔPSI| > 0.15") +
  scale_fill_brewer(palette="Set1")
