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

# Add metadata columns to each dataframe beforehand
dfs_with_metadata <- lapply(dfs, function(x) {
  x$df %>%
    mutate(Group = x$group, Study = x$study)
})

# Combine all dataframes into one using bind_rows() for faster merging
final_df <- bind_rows(dfs_with_metadata) %>%
  select(EVENT, Group, Study, deltapsi) %>%
  filter(abs(deltapsi)>=0.15) %>%
  na.omit()# Adjust this line if other columns should be prioritized


plot_psi_distribution <- ggplot(final_df, aes(x = deltapsi, fill = ifelse(deltapsi < 0, "Skipped", "Included"))) +
  geom_histogram(alpha = 0.8, bins = 100, position = "identity") +
  labs(
    title = "PSI Distribution of all Events",
    subtitle = "All 5 studies were merged",
    x = "ΔPSI (Myotubes-Myoblasts)",
    y = "Density",
    fill = "|ΔPSI| > 0.15",
    caption = paste0("Created by AG on ", Sys.Date())
  ) +
  theme_minimal(base_family = font) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(size = 15),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.3)
  ) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~Group)

