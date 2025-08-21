# Required Packages
library(AzureStor)
library(processx)
library(dplyr)
library(tidyr)
library(glue)
library(stringr)

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Please provide the required input.")
} else {
  sas <- args[1]
  samples_path <- args[2]
  paired <- args[3]
  out_dir <- args[4]
  publish_dir <- args[5]
}

samples <- read.csv(samples_path)

# prep sample to fastq path df
if (paired == "true") {
  cols <- c("sample_id", "read1", "read2")
} else {
  cols <- c("sample_id", "read")  
}
sample_fastq <- data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols))) %>%
  mutate_all(as.character)
sample_fastq_pub <- data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols))) %>%
  mutate_all(as.character)

# loop over each sample and download fastqs
for (i in 1:nrow(samples)) {
  sample_id <- samples$sample_id[i]
  az_url <- samples$az_url[i]

  az_url <- ifelse(grepl("http", az_url), az_url, paste0("https://", az_url))

  AzureStor::call_azcopy("copy",
                          glue("{az_url}/*?{sas}"),
                          glue({"{out_dir}/{sample_id}"}),
                          '--include-regex', ".*\\.(fastq\\.gz|fq\\.gz|fastq|fq)")

  los <- list.files(path = glue("{out_dir}/{sample_id}"), full.names = TRUE) %>%
    normalizePath()

  # rename fastq files with proper sample ID
  for (i in 1:length(los)) {
    exts <- str_extract(los[i], "\\..+")
    dir_path <- dirname(los[i])
    new_path <- glue("{dir_path}/{sample_id}_{i}{exts}")

    file.rename(los[i], new_path)
    los[i] <- new_path
  }

  # prepare output DFs
  los_pub <- str_replace(los, glue("^.*(?=/{out_dir})"), publish_dir)
  row <- data.frame(t(c(sample_id, los)))
  row_pub <- data.frame(t(c(sample_id, los_pub)))
  if (ncol(row) != length(cols)) {
    warning(glue("{sample_id} did not match to the correct number of fastq file(s) (expects {length(cols)-1}, but matched {ncol(row) - 1})"))
    next
  }
  colnames(row) <- cols
  colnames(row_pub) <- cols
  sample_fastq <- bind_rows(sample_fastq, row)
  sample_fastq_pub <- bind_rows(sample_fastq_pub, row_pub)
}

write.csv(sample_fastq, glue("sample_fastq.csv"), row.names = FALSE, quote = FALSE)
write.csv(sample_fastq_pub, glue("sample_fastq_pub.csv"), row.names = FALSE, quote = FALSE)
