nextflow.enable.dsl=2

include { TRIMMOMATIC } from '../../modules/local/trimmomatic'
include { ALIGNMENT as BWA_MEM2_ALIGNMENT } from '../../modules/local/alignment'
include { SAMTOOLS_SORT_INDEX } from '../../modules/local/samtools_sort_index'
include { MARK_DUPLICATES } from '../../modules/local/mark_duplicates'

workflow ALIGNMENT_PROCESSING {

	take:
		fastq_ch
		paired_end
		reference_genome_path
		reference_genome_version
		crop_50
		output_directory
		
	main:

        TRIMMOMATIC(fastq_ch)
		BWA_MEM2_ALIGNMENT(TRIMMOMATIC.out, reference_genome_path.parent) 	// passing parent dir path as input properly stages directory in container 	
		SAMTOOLS_SORT_INDEX(BWA_MEM2_ALIGNMENT.out)
		MARK_DUPLICATES(SAMTOOLS_SORT_INDEX.out)

	emit:
		// channel of [sample_id, bam, bai, txt] per sample
		id_bams_ch = MARK_DUPLICATES.out[0]
		// channel of bams
		bams_ch = MARK_DUPLICATES.out[1]
}
