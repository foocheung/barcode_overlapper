
base_dir <- "/Volumes/chi-1/PROJECTS/TEMP/240716_VH00286_153_AACCFJWHV/RUN"


# extract_compare_barcodes.R

# Load necessary libraries
library(ShortRead)
library(ggplot2)
library(reshape2)

# Function to extract barcodes from a subsample of a fastq file using streaming
extract_barcodes <- function(fastq_path, barcode_length = 16, subsample_size = 1e6) {
  # Create a FastqStreamer object to stream the FASTQ file
  streamer <- FastqStreamer(fastq_path, n = subsample_size)
  
  # Initialize a vector to store barcodes
  barcodes <- character(0)
  
  # Stream the fastq file in chunks
  while (length(fastq_chunk <- yield(streamer)) > 0) {
    # Extract sequences from the chunk
    sequences <- sread(fastq_chunk)
    
    # Extract the barcodes (first `barcode_length` bases of each sequence)
    chunk_barcodes <- subseq(sequences, start = 1, end = barcode_length)
    
    # Convert to character vector and combine with previously extracted barcodes
    barcodes <- c(barcodes, as.character(chunk_barcodes))
    
    # Stop if we have enough barcodes
    if (length(barcodes) >= subsample_size) {
      barcodes <- barcodes[1:subsample_size]
      break
    }
  }
  
  # Remove duplicates
  unique(barcodes)
}

# Function to calculate barcode overlap
calculate_overlap <- function(gex_barcodes, adt_barcodes) {
  length(intersect(gex_barcodes, adt_barcodes))
}

# Function to find all R1 fastq files in specified folders (GEX$number and ADT$number)
find_fastq_files <- function(base_dir, prefix) {
  dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  target_dirs <- dirs[grepl(paste0("^", prefix, "\\d+$"), basename(dirs))]
  fastq_files <- list.files(target_dirs, pattern = ".*_R1_001.fastq.gz$", recursive = TRUE, full.names = TRUE)
  return(fastq_files)
}



# Find all GEX and ADT fastq files
gex_files <- find_fastq_files(base_dir, "GEX")
adt_files <- find_fastq_files(base_dir, "ADT")

# Extract barcodes for all GEX and ADT files
gex_barcodes_list <- lapply(gex_files, extract_barcodes)
adt_barcodes_list <- lapply(adt_files, extract_barcodes)

# Create a matrix to store overlap counts
overlap_matrix <- matrix(0, nrow = length(gex_files), ncol = length(adt_files))
rownames(overlap_matrix) <- basename(dirname(gex_files))
colnames(overlap_matrix) <- basename(dirname(adt_files))

# Calculate overlaps
for (i in seq_along(gex_files)) {
  for (j in seq_along(adt_files)) {
    overlap_matrix[i, j] <- calculate_overlap(gex_barcodes_list[[i]], adt_barcodes_list[[j]])
  }
}

# Convert the matrix to a long format data frame
overlap_df <- melt(overlap_matrix, varnames = c("GEX", "ADT"), value.name = "Overlap")

# Plot the heatmap using ggplot2 with geom_tile
ggplot(overlap_df, aes(x = ADT, y = GEX, fill = Overlap)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Barcode Overlap Between GEX and ADT Samples",
       x = "ADT Sample",
       y = "GEX Sample",
       fill = "Overlap")

# Save the heatmap to a file
ggsave("barcode_overlap_heatmap.png")

