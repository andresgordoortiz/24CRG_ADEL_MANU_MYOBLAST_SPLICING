# Create clust tables
clust_events_tao<-na.omit(data.frame(ID=tao_events$PSI$EVENT, t0=rowMeans(tao_events$PSI[,7:9]), t2=rowMeans(tao_events$PSI[,10:12]),t4=rowMeans(tao_events$PSI[,13:15]),t7=rowMeans(tao_events$PSI[,16:18])))
write_tsv(clust_events_tao, "clust_events_tao.txt")

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

