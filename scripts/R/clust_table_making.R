# Create clust tables gene expression

splicing_factor_genes<-txi$counts[rownames(txi$counts) %in% mouse_sf_filtered,]
colnames(splicing_factor_genes)<-metadata$run_accession
clust_events_tao<-na.omit(data.frame(ID=rownames(splicing_factor_genes), t0=rowMeans(splicing_factor_genes[,1:3]), t2=rowMeans(splicing_factor_genes[,4:6]),t4=rowMeans(splicing_factor_genes[,7:9]),t7=rowMeans(splicing_factor_genes[,10:12])))
write_tsv(clust_events_tao, "clust_events_tao_gene_expression.txt")

clust_events_trapnell<-na.omit(data.frame(ID=trapnell_events$PSI$EVENT, t0=rowMeans(trapnell_events$PSI[,7:9]), t2=rowMeans(trapnell_events$PSI[,10:11]),t4=rowMeans(trapnell_events$PSI[,12:13]),t7=rowMeans(trapnell_events$PSI[,14:16])))
write_tsv(clust_events_trapnell, "clust_events_trapnell.txt")

# Read clust output


# Read the TSV file
# Skip the first row (for descriptive cluster information) and specify the second row as column names
clust_out <- read.delim("Clusters_Objects.tsv", header = FALSE, skip = 1, stringsAsFactors = FALSE)[-1,]

# Set proper column names using the first row of the file
column_names <- read.delim("Clusters_Objects.tsv", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
colnames(clust_out) <- column_names
clust_out[clust_out == ""] <- NA

