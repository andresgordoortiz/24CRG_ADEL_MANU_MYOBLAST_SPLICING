# Create clust tables gene expression



#clust_events_tao<-na.omit(data.frame(ID=rownames(splicing_factor_genes), t0=rowMeans(splicing_factor_genes[,1:3]), t2=rowMeans(splicing_factor_genes[,4:6]),t4=rowMeans(splicing_factor_genes[,7:9]),t7=rowMeans(splicing_factor_genes[,10:12])))
write_tsv(splicing_factor_genes, "clust_events_tao_splicing.txt")

clust_events_trapnell<-na.omit(data.frame(ID=trapnell_events$PSI$EVENT, t0=rowMeans(trapnell_events$PSI[,7:9]), t2=rowMeans(trapnell_events$PSI[,10:11]),t4=rowMeans(trapnell_events$PSI[,12:13]),t7=rowMeans(trapnell_events$PSI[,14:16])))
write_tsv(clust_events_trapnell, "clust_events_trapnell.txt")

splicing_factor_genes<-tao_events$PSI[grep("IN",tao_events$PSI$EVENT),c(2,7:18)]
splicing_factor_genes$ID<-splicing_factor_genes$EVENT
splicing_factor_genes <- splicing_factor_genes[, c("ID",metadata_tao$run_accession)]
clust_events_tao<-na.omit(data.frame(ID=splicing_factor_genes$ID, t0=rowMeans(splicing_factor_genes[,2:4]), t2=rowMeans(splicing_factor_genes[,5:7]),t4=rowMeans(splicing_factor_genes[,8:10]),t7=rowMeans(splicing_factor_genes[,11:13])))
write_tsv(clust_events_tao, "clust_events_splicing_tao_onlyintrons.txt")


splicing_factor_genes<-txi$counts
splicing_factor_genes$ID<-rownames(splicing_factor_genes)
splicing_factor_genes <- splicing_factor_genes[, c("ID",metadata_tao$run_accession)]
clust_events_tao<-na.omit(data.frame(ID=splicing_factor_genes$ID, t0=rowMeans(splicing_factor_genes[,2:4]), t2=rowMeans(splicing_factor_genes[,5:7]),t4=rowMeans(splicing_factor_genes[,8:10]),t7=rowMeans(splicing_factor_genes[,11:13])))
write_tsv(clust_events_tao, "clust_events_splicing_tao.txt")
# Read clust output


# Read the TSV file
# Skip the first row (for descriptive cluster information) and specify the second row as column names
clust_out <- read.delim("Clusters_Objects.tsv", header = FALSE, skip = 1, stringsAsFactors = FALSE)[-1,]

# Set proper column names using the first row of the file
column_names <- read.delim("Clusters_Objects.tsv", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
colnames(clust_out) <- column_names
clust_out[clust_out == ""] <- NA

clustered_exons<-read_tsv("results/clust_results_tao_exons/Clusters_Objects.tsv")[-1,]