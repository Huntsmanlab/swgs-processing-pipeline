nextflow.enable.dsl=2

include { AZURE_DOWNLOAD } from '../../modules/local/azure_download'

workflow COLLECT_FASTQS {

	take:
		from_azure
		fastq_csv
		paired_end


	main:

        if (from_azure) {

            AZURE_DOWNLOAD(fastq_csv, paired_end)
			AZURE_DOWNLOAD.out[1].view()
            fastq_ch = AZURE_DOWNLOAD.out[1] \
                | splitCsv(skip: 1) \
                | map { row -> [row[0], row[1..-1]] }

        } else {

            fastq_ch = Channel.fromPath(fastq_csv, checkIfExists: true) \
                | splitCsv(skip: 1) \
                | map { row -> [row[0], row[1..-1]] }
        }

	emit:
		fastq_ch
}
