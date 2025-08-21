#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(future)
  library(CGHcall)
  library(data.table)
})

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Please provide the required input.")
} else {
  nthreads <- args[1]
  type <- args[2]
  dir_path <- args[3]
  genome <- args[4]
  out_path <- args[5]
}

binsize <- as.numeric(paste0(stringr::str_remove(basename(out_path), 'kb'), "000"))

future::plan("multisession", workers=as.integer(nthreads))

write_plot_acn <- function(sample, path, obj) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  plot <- utanos::CNSegmentsPlot(obj, sample = sample, def_point_colour = 'grey',
                                 copy_number_breaks = c(2,4,6,8,10,12),
                                 max_copy_number = 15) + ggplot2::ggtitle(sample)
  ggplot2::ggsave(paste0(path, "/", sample, ".png"), plot, width = 10, height = 6, scale = 2)
}

# Convert WX files to QDNAseq object
if (type == "wx") {
  cnobj <- utanos::BuildWxQdnaObject(dir_path, genome_used = genome, 
                                       bin_size = binsize, Xincluded = TRUE)
  
  # find solutions
  solutions <- utanos::FindRascalSolutions(cnobj)
  write.csv(solutions, glue::glue("{out_path}/{type}_{binsize/1000}kb_rascal_solutions.csv"), row.names = FALSE)
  
  # scaling to ACN
  results <- utanos::CalculateACNs(cnobj, "mad", solutions, return_sols = TRUE, return_S4 = TRUE)
  write.csv(results[['rascal_solutions']], glue::glue("{out_path}/{type}_{binsize/1000}kb_chosen_solns.csv"), row.names = FALSE)
  acn_obj <- results[["acn.obj"]]
  saveRDS(acn_obj, glue::glue("{out_path}/{type}_{binsize/1000}kb_aCN_segments.rds"))

  # plots 
  sample_names <- sampleNames(acn_obj)
  lapply(sample_names, write_plot_acn, path = glue::glue("{out_path}/acn_plots"), obj = acn_obj)

} else if (type == "qdnaseq") {
  for (i in c("noY","noXY")) {
    cnobj <- readRDS(glue::glue("{dir_path}/{basename(dir_path)}_{i}_copyNumbersSegmented.rds"))
    
    # find solutions
    solutions <- utanos::FindRascalSolutions(cnobj)
    write.csv(solutions, glue::glue("{out_path}/{type}_{binsize/1000}kb_rascal_solutions_{i}.csv"), row.names = FALSE)
    
    # scaling to ACN
    results <- utanos::CalculateACNs(cnobj, "mad", solutions, return_sols = TRUE, return_S4 = TRUE)
    write.csv(results[['rascal_solutions']], glue::glue("{out_path}/{type}_{binsize/1000}kb_chosen_solns_{i}.csv"), row.names = FALSE)
    acn_obj <- results[["acn.obj"]]
    saveRDS(acn_obj, glue::glue("{out_path}/{type}_{binsize/1000}kb_aCN_segments_{i}.rds"))
    
    # plots
    sample_names <- sampleNames(acn_obj)
    lapply(sample_names, write_plot_acn, path = glue::glue("{out_path}/acn_plots/{i}"), obj = acn_obj)  
  }
}
