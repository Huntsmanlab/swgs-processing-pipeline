nextflow.enable.dsl=2

include { CN_ICHOR } from '../../modules/local/cn_ichor'
include { ICHOR_PROCESS } from '../../modules/local/ichor_process'
include { ICHOR_COMBINE } from '../../modules/local/ichor_combine'

workflow ICHOR {

	take:
		id_bams_ch
		bins_ch
		paired_end

	main:

        if (paired_end) {

			bam_type = Channel.from("pe")
            bam_ch = id_bams_ch.combine(bam_type)

        } else {

            bam_type = Channel.from("se")
            bam_ch = id_bams_ch.combine(bam_type)
        }		

		ichor_input = bam_ch.combine(bins_ch)
		bam_ch_filtered = ichor_input.filter { val -> val[5] in [10, 50, 500, 1000] }
		
		CN_ICHOR(bam_ch_filtered)
		ICHOR_PROCESS(CN_ICHOR.out.results)	
		
		ichor_collect = ICHOR_PROCESS.out.groupTuple(by: [0,1])
		ICHOR_COMBINE(ichor_collect)
}
