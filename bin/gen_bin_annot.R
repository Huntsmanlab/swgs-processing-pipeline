#!/usr/bin/env Rscript

# Required packages
library(QDNAseq)
library(Biobase)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Please provide the required input.")
} else {
  nthreads <- args[1]
  bin_size <- as.numeric(args[2])
  genome <- args[3]
  mappability <- args[4]
  blacklist <- args[5]
  bigwig_average <- args[6]
  bam_path <- args[7]
  out_path <- args[8]
  paired <- args[9]
}

# Setup for parallelization
options(future.globals.maxSize = 4 * 1024^3) # Set max global size to 4 GB
plan(strategy = "multicore", workers = as.numeric(nthreads) - 1) # Workers = cores - 1; adjust as appropriate for your Slurm job

if (genome == 'mm10') {
  bsgenome = BSgenome.Mmusculus.UCSC.mm10
  chrs <- 1:19
} else if (genome == 'hg19') {
  bsgenome = BSgenome.Hsapiens.UCSC.hg19
  chrs <- 1:22
} else if (genome == 'hg38') {
  bsgenome = BSgenome.Hsapiens.UCSC.hg38
  chrs <- 1:22
} else {
  stop("Invalid genome.")
}

bins <- createBins(bsgenome = bsgenome, binSize = bin_size)

# Get the mappability per bin
mappaCol <- calculateMappability(bins = bins, bigWigFile = mappability, bigWigAverageOverBed = bigwig_average)

if (sum(mappaCol) == 0) {
  mappaCol <- calculateMappability(bins = bins, bigWigFile = mappability, bigWigAverageOverBed = bigwig_average, chrPrefix = "")
}

bins$mappability <- mappaCol

# Calculate bin overlap with the blacklisted regions
bins$blacklist <- calculateBlacklist(bins = bins, bedFiles = blacklist)

# Get the control BAM files and calculate the LOESS residuals per bin; exclude the mitochondrial genome and bins with all N bases
control_samples <- binReadCounts(bins = bins, path = bam_path, ext = 'se.*bam', pairedEnds = FALSE)
control_samples <- applyFilters(object = control_samples, residual = FALSE, blacklist = FALSE, mappability = FALSE, bases = FALSE, chromosomes = c("MT", "Y"))

bins$residual <- iterateResiduals(object = control_samples) # Populates the residuals for each bin with the median residual of the control samples
bins$use <- bins$chromosome %in% as.character(chrs) & bins$bases > 0 # By default, use all bins that are autosomal and have at least one characterized base

# Save our annotated bins for later use
saveRDS(bins, file = paste0(out_path, "/", genome, "_bins_", bin_size, "kbp_SE150.rds"))

if (paired == "true") {
  # Get the control BAM files and calculate the LOESS residuals per bin; exclude the mitochondrial genome and bins with all N bases
  control_samples <- binReadCounts(bins = bins, path = bam_path, ext = 'pe.*bam', pairedEnds = TRUE)
  control_samples <- applyFilters(object = control_samples, residual = FALSE, blacklist = FALSE, mappability = FALSE, bases = FALSE, chromosomes = c("MT", "Y"))

  bins$residual <- iterateResiduals(object = control_samples) # Populates the residuals for each bin with the median residual of the control samples
  bins$use <- bins$chromosome %in% as.character(chrs) & bins$bases > 0 # By default, use all bins that are autosomal and have at least one characterized base

  # Save our annotated bins for later use
  saveRDS(bins, file = paste0(out_path, "/", genome, "_bins_", bin_size, "kbp_PE150.rds"))
}
