nextflow.enable.dsl=2

include { AZURE_DOWNLOAD } from '../../modules/local/azure_download'

workflow COLLECT_BAMS {

	take:
		fastq_csv
		from_azure
		paired_end
		
	main:

        bams_ch = Channel.fromPath(fastq_csv, checkIfExists: true) \
            | splitCsv(skip: 1) \
            | map { row -> row[1] }

        id_bams_ch = Channel.fromPath(fastq_csv, checkIfExists: true) \
            | splitCsv(skip: 1) \
            | map { row -> [row[0], row[1], row[2], row[3]] }

	emit:	
		id_bams_ch
		bams_ch
}
