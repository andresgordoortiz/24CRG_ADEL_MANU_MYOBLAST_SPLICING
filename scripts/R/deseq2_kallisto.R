#!/usr/bin/env Rscript

# -------------------------------
# Load required libraries
# -------------------------------
if (!requireNamespace("tximport", quietly = TRUE))
  install.packages("tximport")
if (!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("biomaRt")
if (!requireNamespace("readr", quietly = TRUE))
  install.packages("readr")
if (!requireNamespace("DESeq2", quietly = TRUE))
  install.packages("DESeq2")

library(tximport)
library(biomaRt)
library(readr)
library(DESeq2)


# -------------------------------
# Define file paths and read sample metadata
# -------------------------------
# This sample table should have columns: sample and condition
sampleTableFile <- arrange(read.csv(paste0(getwd(),"/metadata/tao_metadata.csv"), sep = ","), experiment_title)

sampleTableFile$experiment_title <- factor(sampleTableFile$experiment_title)

# Connect to Ensembl biomart for mouse genes
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve transcript-to-gene mapping (ENSMUST to ENSMUSG)
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name"),
                 mart = mart)

# Rename columns to match tximport format
colnames(tx2gene) <- c("TXNAME", "GENEID")

# -------------------------------
# Define file paths for kallisto output files
# -------------------------------
# Assumes that each sample's kallisto results are in "kallisto_outputs/<sample>/abundance.tsv"
files <- file.path("notebooks/final/Kallisto", sampleTableFile$run_accession, "abundance.tsv")
names(files) <- sampleTableFile$run_accession

# Check that all files exist
if (!all(file.exists(files))) {
  missingFiles <- files[!file.exists(files)]
  stop("The following kallisto files were not found:\n", paste(missingFiles, collapse = "\n"))
}

# -------------------------------
# Import kallisto quantification data using tximport
# -------------------------------
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# -------------------------------
# Create DESeq2 dataset from the tximport object
# -------------------------------
dds <- DESeqDataSetFromTximport(txi, colData = sampleTableFile, design = ~ experiment_title)

# Optionally, pre-filter low-count genes (e.g., keep rows with at least 10 counts across samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Extract the fourth quantile of expression
# Calculate the mean TPM for each gene across samples
meanTPM <- rowMeans(txi$abundance)

# Plot the distribution of mean TPM values
hist(meanTPM, breaks = 50, main = "Distribution of Mean TPM", xlab = "Mean TPM")

# Determine the threshold for the top 25% highest expressed genes
threshold <- quantile(meanTPM, probs = 0.75)

# Select genes with mean TPM above this threshold
highExpGenes <- names(meanTPM)[meanTPM >= threshold]

saveRDS(highExpGenes, file = "highExpGenes.rds")

# Determine the threshold for the top 25% lowest
threshold_low <- quantile(meanTPM, probs = 0.25)

lowExpGenes <- names(meanTPM)[meanTPM <= threshold_low]

saveRDS(lowExpGenes, file = "lowExpGenes.rds")

# -------------------------------
# Run the differential expression analysis with DESeq2
# -------------------------------
dds <- DESeq(dds)

res_d4_vs_d0 <- results(dds, contrast = c("experiment_title", "4_days","0_days"))

df <- as.data.frame(res_d4_vs_d0)
df <- cbind(Gene = rownames(df), df)  # Add row names as a new column
write.csv(df, file = "DESeq2_d4_vs_d0.csv", row.names = FALSE, sep = ",")  

res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Save the full results table to a CSV file
write.csv(as.data.frame(resOrdered), file = "DESeq2_results.csv", row.names = TRUE)

# -------------------------------
# Apply a variance stabilizing transformation (VST)
# -------------------------------
# VST is preferred for visualizations because it makes the data more homoscedastic
vsd <- vst(dds, blind = FALSE)

# -------------------------------
# PCA Plot
# -------------------------------
pcaData <- plotPCA(vsd, intgroup = "experiment_title", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = experiment_title)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal()
print(p)

# -------------------------------
# Heatmap of Sample Distances
# -------------------------------
# Calculate sample-to-sample distances using the VST-transformed data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap",
         filename = "sample_distance_heatmap.png")
