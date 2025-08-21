nextflow.enable.dsl=2

// validate parameter inputs
def assert_required_param(param, param_name){
    if(param == null){
        exit 1, param_name +' not specified. Please provide ${param_name}!'
    } else {
		println "${param_name} = ${param}"
	}
}
assert_required_param(params.input_csv, 'input_csv')
assert_required_param(params.from_azure, 'from_azure')
assert_required_param(params.paired_end, 'paired_end')
assert_required_param(params.reference_genome_path, 'reference_genome_path')
assert_required_param(params.reference_genome_version, 'reference_genome_version')
assert_required_param(params.crop_50, 'crop_50')
assert_required_param(params.input_directory, 'input_directory')
assert_required_param(params.output_directory, 'output_directory')
assert_required_param(params.input_bam, 'input_bam')
assert_required_param(params.run_qdnaseq, 'run_qdnaseq')
assert_required_param(params.bin_sizes, 'bin_sizes')
assert_required_param(params.create_bin_annotations, 'create_bin_annotations')
assert_required_param(params.run_ichor, 'run_ichor')
assert_required_param(params.run_wisecondorx, 'run_wisecondorx')

input_csv = file(params.input_csv)
from_azure = params.from_azure
paired_end = params.paired_end
reference_genome_path = file(params.reference_genome_path)
reference_genome_version = params.reference_genome_version
crop_50 = params.crop_50
input_directory = params.input_directory
output_directory = params.output_directory
input_bam = params.input_bam
bins_ch = Channel.from(params.bin_sizes)

run_qdnaseq = params.run_qdnaseq
create_bin_annotations = params.create_bin_annotations
run_wisecondorx = params.run_wisecondorx
run_ichor = params.run_ichor
bin_annotations = params.bin_annotations
wcx_ref = params.wcx_ref
wcx_create_ref = params.wcx_create_ref

// import subworkflows
include { COLLECT_FASTQS } from '../subworkflows/local/collect_fastqs'
include { COLLECT_BAMS } from '../subworkflows/local/collect_bams'
include { PRE_ALIGNMENT_QC } from '../subworkflows/local/pre_alignment_qc'
include { ALIGNMENT_PROCESSING } from '../subworkflows/local/alignment_processing'
include { POST_ALIGNMENT_QC } from '../subworkflows/local/post_alignment_qc'
include { QDNASEQ as CN_QDNASEQ } from '../subworkflows/local/cn_qdnaseq' 
include { ICHOR as CN_ICHOR } from '../subworkflows/local/cn_ichor'
include { WISECONDORX as CN_WISECONDORX } from '../subworkflows/local/cn_wisecondorx'

// run main workflow
workflow SWGS_QC_PIPELINE {

	if (input_bam) {
		
		COLLECT_BAMS(
			input_csv,
			from_azure,
			paired_end
		)

		id_bams_ch = COLLECT_BAMS.out.id_bams_ch
		bams_ch = COLLECT_BAMS.out.bams_ch

	} else {

		COLLECT_FASTQS(
        	from_azure,
        	input_csv,
        	paired_end
		)
	
		fastq_ch = COLLECT_FASTQS.out
		PRE_ALIGNMENT_QC(fastq_ch)

		ALIGNMENT_PROCESSING(
        	fastq_ch,
        	paired_end,
        	reference_genome_path,
        	reference_genome_version,
        	crop_50,
        	output_directory
		)

        id_bams_ch = ALIGNMENT_PROCESSING.out.id_bams_ch
        bams_ch = ALIGNMENT_PROCESSING.out.bams_ch
	}

	POST_ALIGNMENT_QC(id_bams_ch)
	
	if (run_qdnaseq || create_bin_annotations) {
		CN_QDNASEQ(
			bams_ch,
			bins_ch,
			create_bin_annotations,
			bin_annotations,
			paired_end,
			reference_genome_path,
			reference_genome_version,
			run_qdnaseq
		)
	}

	if (run_ichor) {
		CN_ICHOR(
			id_bams_ch,
			bins_ch,
			paired_end
		)
	}

	if (run_wisecondorx) {
		CN_WISECONDORX(
			bams_ch,
			bins_ch,
			paired_end,
			wcx_ref,
			wcx_create_ref
		)		
	}
}
