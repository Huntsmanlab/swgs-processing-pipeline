process CN_QDNASEQ_RCN {
    publishDir params.output_directory, mode:'copy'
    tag "Completing relative copy number analysis via QDNAseq"
    container'docker://huntsmanlab/qdnaseq:latest' 
    cpus 16
    memory '32 G'
    time '12h'

    input:
    tuple val(bin_size), val(bam_type), val(bams)
	path(binannos_dir)

    output:
    tuple val(bin_size), val(bam_type)
	path("relative_cn/qdnaseq/${bam_type}/${bin_size}kb")

    script:
    """
    mkdir -p "relative_cn/qdnaseq/${bam_type}/${bin_size}kb"
	printf '%s\n' "${bams.join('\n')}" > bamfileslist.txt
    shopt -s nocaseglob extglob
   	runQDNAseq.R ${bin_size}kb ${task.cpus} relative_cn/qdnaseq/${bam_type}/${bin_size}kb bamfileslist.txt ${binannos_dir}/*@(${bin_size}kb|${bam_type})*@(${bin_size}kb|${bam_type})*.rds
    rm bamfileslist.txt
    """
}
