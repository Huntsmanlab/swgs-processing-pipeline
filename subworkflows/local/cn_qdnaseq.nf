nextflow.enable.dsl=2

include { QDNA_BINS as GEN_QDNA_ANNOTS } from '../../modules/local/qdna_bin_annot'
include { CN_QDNASEQ_RCN } from '../../modules/local/cn_qdnaseq_rcn'
include { CN_QDNASEQ_ACN } from '../../modules/local/cn_qdnaseq_acn'

workflow QDNASEQ {

	take:
		bam_ch
		bins_ch
		create_bin_annotations
		bin_annotations
		paired_end
		reference_genome_path
		reference_genome_version
		run_qdnaseq

	main:
		
		ref_genome_ch = channel.from(reference_genome_version)

		if (paired_end) {
			bam_ls = bam_ch.map { bam -> ["pe", bam] }

		} else {
            bam_ls = bam_ch.map { bam -> ["se", bam] }
		}

		// generate QDNAseq bin annotations
		if (create_bin_annotations) {
			bin_annotations = GEN_QDNA_ANNOTS(bins_ch, reference_genome_path, paired_end)
			qdnaseq_ch = bins_ch.combine(bam_ls).flatten()

		} else {
			qdnaseq_ch = bins_ch.combine(bam_ls) | groupTuple(by: [0,1])
		}

		// CN analysis via QDNAseq
		if (run_qdnaseq) {
		
			rcn_ch = CN_QDNASEQ_RCN(qdnaseq_ch, bin_annotations)
			qdna_ch = rcn_ch[0].combine(ref_genome_ch)
			CN_QDNASEQ_ACN(qdna_ch, rcn_ch[1])
		}
}
