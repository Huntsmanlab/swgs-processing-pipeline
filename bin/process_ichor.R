#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(dplyr)
library(data.table)
})

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
cat("Number of arguments received:", length(args), "\n")

# Check if arguments are provided
if (length(args) < 5) {
  stop("Please provide the required arguments")
} else {
  sample_id <- args[1]
  bins_file <- args[2]
  segs_file <- args[3]
  out_dir <- args[4]
  genome <- args[5]
}

# read in relative CN per bins and segments
bins_df <- read.delim(file = bins_file, sep = "\t")
segs_df <- read.delim(file = segs_file, sep = "\t")

bins_df_clean <- bins_df %>%
  dplyr::rename("ratio" = "log2_TNratio_corrected", "chromosome" = "chr")

segs_df_clean <- segs_df %>%
  dplyr::rename("segVal" = "median", "sample_id" = "sample", "chromosome" = "chr") %>%
  dplyr::select(sample_id, chromosome, start, end, segVal)

# determine bin size
bin_size = bins_df_clean$end[1] - bins_df_clean$start[1] + 1

# transform segment CN table to bin CN table
segs_df_clean_binned <- utanos::SegmentsToCopyNumber(segs = segs_df_clean, bin_size = bin_size, genome = genome, Xincluded = TRUE) %>%
  dplyr::rename("ratio_median" = "segmented")

# combine tables
bins_segs_df <- inner_join(x = bins_df_clean, y = segs_df_clean_binned, by = join_by(chromosome, start, end)) %>%
  dplyr::select(sample_id, everything()) %>%
  dplyr::select(!end)

# create output file
write.table(x = bins_segs_df, paste0(out_dir, sample_id, "_rCN.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
