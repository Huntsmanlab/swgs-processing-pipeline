nextflow.enable.dsl=2

multiqc_config = file("$projectDir/config/multiqc_config.yaml")

include { SAMTOOLS_COV_FLAGSTAT } from '../../modules/local/samtools_cov_flagstat'
include { FASTQC_POSTALIGN } from '../../modules/local/fastqc_postalign'
include { MULTIQC_POSTALIGN } from '../../modules/local/multiqc_postalign'

workflow POST_ALIGNMENT_QC {

	take:
		id_bams_ch

	main:
	
		SAMTOOLS_COV_FLAGSTAT(id_bams_ch)
		FASTQC_POSTALIGN(id_bams_ch)

		path_ch = id_bams_ch.map { val -> [val[1], val[2], val[3]] }
		multiqc_logs = path_ch.mix(SAMTOOLS_COV_FLAGSTAT.out.log1, SAMTOOLS_COV_FLAGSTAT.out.log2, FASTQC_POSTALIGN.out.log).collect()

        multiqc_config_ch = Channel.fromPath(multiqc_config, checkIfExists: true)
        MULTIQC_POSTALIGN(multiqc_logs, multiqc_config_ch)
}
