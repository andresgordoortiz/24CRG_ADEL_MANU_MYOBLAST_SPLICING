# Load necessary Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("BSgenome", "BSgenome.Mmusculus.UCSC.mm10", "GenomicRanges", "Biostrings", "rtracklayer", "ggplot2", "Gviz"))

# Load libraries
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)
library(ggplot2)
library(Gviz)

# Define the genome object
genome <- BSgenome.Mmusculus.UCSC.mm10

# Load a set of exon coordinates (example: GTF file or predefined GRanges)
# Replace this with your own exon data or load from a file
# Example: exons <- rtracklayer::import("path_to_exons.gtf")

mef2d_MmuEX0028567_MmuEX0028566 <- GRanges(seqnames = c("chr3"), 
                              ranges = IRanges(start = 88156322, end = 88159332),
                              strand = c("+"))


# Extract sequences of flanking regions
event_seqs <- getSeq(genome, mef2d_MmuEX0028567_MmuEX0028566)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)

rapgef1_MmuEX0038772_MmuEX0038771 <- GRanges(seqnames = c("chr2","chr2"), 
                                             ranges = IRanges(start = c(29710153, 29707198), end = c(29722309, 29722309)),
                                             strand = c("+", "+"))

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, rapgef1_MmuEX0038772_MmuEX0038771)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)


pdlim7_MmuEX0034296_MmuEX0034298 <- GRanges(seqnames = c("chr13","chr13"), 
                                             ranges = IRanges(start = c(55506306, 55507488), end = c(55507366, 55508326)),
                                             strand = c("-", "-"))

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, pdlim7_MmuEX0034296_MmuEX0034298)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)


fermt2_MmuEX0019078 <- GRanges(seqnames = "chr14", 
                               ranges = IRanges(start = 45475888, end = 45482566),
                               strand =  "-")

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, fermt2_MmuEX0019078)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)


arlbp2_MmuINT0018915 <- GRanges(seqnames = "chr8", 
                               ranges = IRanges(start = 94666600, end = 94667658),
                               strand =  "+")

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, arlbp2_MmuINT0018915)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)

araf_MmuINT1004569 <- GRanges(seqnames = "chrX", 
                                ranges = IRanges(start = 20851930, end = 20853459),
                                strand =  "+")

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, araf_MmuINT1004569)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)

capzb_MmuEX0009114 <- GRanges(seqnames = "chr4", 
                              ranges = IRanges(start = 139287756, end = 139291818),
                              strand =  "+")

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, capzb_MmuEX0009114)

# Define the motif of interest
motif <- DNAString("TGCATG")  # Use RNAString if needed

vcountPattern(motif, event_seqs)
vmatchPattern(motif, event_seqs)

# Define the genomic regions
fxr1_regions <- GRanges(
  seqnames = c("chr3", "chr3", "chr3"),
  ranges = IRanges(
    start = c(34064119, 34064119, 34058024),
    end = c(34070322, 34070322, 34058231)
  ),
  strand = c("+", "+", "+")
)

# Extract sequences of flanking regions
event_seqs <- getSeq(genome, rapgef1_MmuEX0038772_MmuEX0038771)

# Define the motif of interest
motif <- DNAString("TGCATG")

# Find motif matches and convert to GRanges
motif_matches_list <- lapply(seq_along(event_seqs), function(i) {
  match_data <- as.data.frame(matchPattern(motif, event_seqs[[i]]))
  if (nrow(match_data) > 0) {
    GRanges(
      seqnames = seqnames(rapgef1_MmuEX0038772_MmuEX0038771)[i],
      ranges = IRanges(
        start = start(rapgef1_MmuEX0038772_MmuEX0038771)[i] + match_data$start - 1,
        end = start(rapgef1_MmuEX0038772_MmuEX0038771)[i] + match_data$end - 1
      ),
      strand = strand(rapgef1_MmuEX0038772_MmuEX0038771)[i]
    )
  } else {
    NULL
  }
})

# Combine into a single GRanges object, removing NULLs
motif_matches <- do.call(c, motif_matches_list[!sapply(motif_matches_list, is.null)])
if (length(motif_matches) == 0) {
  stop("No motif matches found in the provided sequences.")
}

# Create genomic tracks for visualization
itrack <- IdeogramTrack(genome = "mm10", chromosome = "chr2")
gtrack <- GenomeAxisTrack()
atrack <- AnnotationTrack(
  range = rapgef1_MmuEX0038772_MmuEX0038771,
  name = "Exons/Introns",
  genome = "mm10",
  feature = c("Exon", "Exon")
)
mtrack <- AnnotationTrack(
  range = motif_matches,
  name = "Motif Matches",
  genome = "mm10",
  col = "red",
  fill = "red"
)

# Plot tracks
plotTracks(
  list(itrack, gtrack, atrack, mtrack),
  from = 34058000,
  to = 34071000,
  main = "Gene Locus with Exons and Motif Hits"
)
