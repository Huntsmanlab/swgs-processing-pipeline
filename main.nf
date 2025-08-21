#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SWGS_QC_PIPELINE } from './workflows/swgs_qc'

log.info """\

 sWGS PROCESSING PIPELINE  v 3.0
================================
output directory = ${params.output_directory}
reference genome = ${params.reference_genome_version}
"""

workflow {
	SWGS_QC_PIPELINE()
}
