nextflow.enable.dsl=2

process QDNA_BINS {
	publishDir params.output_directory, mode:'copy'

    input:
    val(bin_size)
	path(ref_genome)
	val(paired_end)

    output:
    path("qdna_bin_annots")

    script: 
    """
    mkdir -p qdna_bin_annots
    WORK_DIR=\$(pwd)
    Rscript gen_bin_annot.R ${task.cpus} ${bin_size} ${ref_genome} ${params.qd_mappability} ${params.qd_blacklist} \$WORK_DIR/${qd_bwavgbed} ${params.qd_nbams} qdna_bin_annots ${paired_end}
    """
}
