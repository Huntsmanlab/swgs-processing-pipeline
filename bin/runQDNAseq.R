#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  # library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(future)
  library(CGHcall)
  library(glue)
})

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments (binsize and output path) are provided
if (length(args) < 1) {
    stop("Please provide the required input.")
} else {
    binsize <- args[1]
    nthreads <- as.integer(args[2])
    output_path <- args[3]
    bamfileslist <- args[4]
    bin_annos <- args[5]
}

future::plan("multisession", workers=nthreads)
# set max total size of defined global variables to 500 GiB
options(future.globals.maxSize = 500 * 1024^3)

# Read-in bams list
# assume bams have been saved to output_bams folder 
bamfiles <- readLines(bamfileslist)
# bam_file_list <- list.files('output_bams')
# bam_file_list <- bam_file_list[grep(".bam$", bam_file_list, ignore.case = TRUE)]

# Get sample names
extract_sample_names <- function(paths) {
  file_names <- gsub('.*/', '', paths)
  pieces_before_period <- gsub('\\..*', '', file_names)
  return(pieces_before_period)
}
sample_names <- extract_sample_names(bamfiles)

# Declare paths
# bin_annos <- paste0('1000genomes_PE150_bin_annotations/hg19_bins_150bp_', binsize, '_SE.rds')

# Load Bin Annotations
print(getwd())
print(list.files())
bins <- readRDS(file = bin_annos)

# Read in bams
readCounts <- binReadCounts(bins, 
                            bamfiles = bamfiles, 
                            bamnames = sample_names)

#### With X-chromosome
#### QDNAseq processing and CN calling
readCountsFiltered <- applyFilters(readCounts, chromosomes = c("Y"), 
                                    residual = TRUE, blacklist = TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

# QDNAseq Segmentation
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

# Call Gains and Losses
glBins <- callBins(copyNumbersSegmented, ncpus = nthreads, method = "cutoff", cutoffs = c(-0.1, 0.1))
cgh_obj <- makeCgh(glBins)

# Fix sample names
copyNumbersSegmented@phenoData@data[["name"]] <- word(copyNumbersSegmented@phenoData@data[["name"]], 1, sep = "\\.")

# Save results
saveRDS(copyNumbersSegmented, file = file.path(output_path, paste0(binsize, "_noY_copyNumbersSegmented.rds")))
saveRDS(cgh_obj, file = file.path(output_path, paste0(binsize, "_noY_gl_rCN.rds")))
saveRDS(glBins, file = file.path(output_path, paste0(binsize, "_noY_bins_rCN.rds")))

# Generate plots
write_plot_rcn <- function(sample, path, obj) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  png(paste0(path, "/", sample, ".png"), width = 2000, height = 1200, pointsize = 16)
  QDNAseq::plot(obj[,sample])
  dev.off()
}

sample_names <- sampleNames(copyNumbersSegmented)
print(glue("generating plots for {binsize} noY"))
lapply(sample_names, write_plot_rcn, path = glue("{output_path}/rcn_plots/noY"), obj = copyNumbersSegmented)


#### Without X-chromosome
#### QDNAseq processing and CN calling
readCountsFiltered <- applyFilters(readCounts, chromosomes = c("X", "Y"), 
                                    residual = TRUE, blacklist = TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

# QDNAseq Segmentation
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

# Call Gains and Losses
glBins <- callBins(copyNumbersSegmented, ncpus = nthreads, method = "cutoff", cutoffs = c(-0.1, 0.1))
cgh_obj <- makeCgh(glBins)

# Fix sample names
copyNumbersSegmented@phenoData@data[["name"]] <- word(copyNumbersSegmented@phenoData@data[["name"]], 1, sep = "\\.")

# Save results
saveRDS(copyNumbersSegmented, file = file.path(output_path, paste0(binsize, "_noXY_copyNumbersSegmented.rds")))
saveRDS(cgh_obj, file = file.path(output_path, paste0(binsize, "_noXY_gl_rCN.rds")))
saveRDS(glBins, file = file.path(output_path, paste0(binsize, "_noXY_bins_rCN.rds")))

# Generate plots
print(glue("generating plots for {binsize} noXY"))
lapply(sample_names, write_plot_rcn, path = glue("{output_path}/rcn_plots/noXY"), obj = copyNumbersSegmented)
