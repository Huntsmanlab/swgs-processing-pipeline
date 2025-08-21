process CN_QDNASEQ_ACN {
    publishDir params.output_directory, mode:'copy'
    tag "QDNAseq - RCN to ACN"
    container'docker://huntsmanlab/utanos:latest' 
    cpus 8
    memory '16 G'
    time '12h'

    input:
	tuple val(bin_size), val(bam_type), val(reference_genome_version)
	path(rcn_dir)

	output:
	path("absolute_cn/qdnaseq/${bam_type}/${bin_size}kb")

	script:
	"""
	mkdir -p "absolute_cn/qdnaseq/${bam_type}/${bin_size}kb"
	runACN.R ${task.cpus} qdnaseq "${rcn_dir}" "${reference_genome_version}" "absolute_cn/qdnaseq/${bam_type}/${bin_size}kb"
	"""
}
