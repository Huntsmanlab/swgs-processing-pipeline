nextflow.enable.dsl=2

include { WCX_REF_CREATE } from '../../modules/local/wcx_ref_create'
include { CN_WISECONDORX_RCN } from '../../modules/local/cn_wisecondorx_rcn'
include { CN_WISECONDORX_ACN } from '../../modules/local/cn_wisecondorx_acn'

workflow WISECONDORX {

    take:
        bam_ch
        bins_ch
        paired_end
		wcx_ref
		wcx_create_ref

    main:

        if (paired_end) {
            bam_ls = bam_ch.map { bam -> ["pe", bam] }.groupTuple()

        } else {
            bam_ls = bam_ch.map { bam -> ["se", bam] }.groupTuple()
        }

        // generate wisecondorx ref from normals
        if (wcx_create_ref) {
		
			wx_ref_ch = WCX_REF_CREATE(params.wcx_normals, bins_ch, params.wcx_norm_bams) 
			wx_input_ch = wx_ref_ch.combine(bam_ls)
			CN_WISECONDORX_RCN(wx_input_ch)
			CN_WISECONDORX_ACN(CN_WISECONDORX_RCN.out)

        } else {

			wx_input_ch = bins_ch.combine(Channel.of(wcx_ref)).combine(bam_ls)
            CN_WISECONDORX_RCN(wx_input_ch, params.output_directory)
            CN_WISECONDORX_ACN(CN_WISECONDORX_RCN.out)
        }
}
