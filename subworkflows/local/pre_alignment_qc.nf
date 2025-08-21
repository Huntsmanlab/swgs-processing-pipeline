nextflow.enable.dsl=2

include { FASTQC_PREALIGN } from '../../modules/local/fastqc_prealign'
include { MULTIQC_PREALIGN } from '../../modules/local/multiqc_prealign'

workflow PRE_ALIGNMENT_QC {

	take:
		fastq_ch

	main:

        FASTQC_PREALIGN(fastq_ch)
        MULTIQC_PREALIGN(FASTQC_PREALIGN.out.logs.collect())
}
