#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SWGS_QC_PIPELINE } from './workflows/swgs_qc'

workflow {
	log.info """\
		sWGS PROCESSING PIPELINE
		================================
		output directory = ${params.output_directory}
		reference genome = ${params.reference_genome_version}
		"""
		.stripIdent()

	SWGS_QC_PIPELINE()
}
