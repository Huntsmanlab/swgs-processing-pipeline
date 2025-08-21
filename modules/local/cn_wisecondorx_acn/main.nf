process CN_WISECONDORX_ACN {
    publishDir params.output_directory, mode:'copy'
    tag "WisecondorX - RCN to ACN w/ bin size ${binsize}kb"
    container'docker://huntsmanlab/utanos:latest' 
    cpus 8
    memory '16 G'
    time '12h'

    input:
    path(rcn_dir)
    val(binsize)
    val(bam_type)

    output:
    path("absolute_cn/wisecondorx/${bam_type}/${binsize}kb")

    script:
    """
	echo ${rcn_dir}
    mkdir wx_files
	cp ${params.output_directory}/relative_cn/wisecondorx/${bam_type}/${binsize}kb/*/*.txt wx_files
    cp ${params.output_directory}/relative_cn/wisecondorx/${bam_type}/${binsize}kb/*/*.bed wx_files
    mkdir -p "absolute_cn/wisecondorx/${bam_type}/${binsize}kb"
    runACN.R ${task.cpus} wx wx_files ${params.reference_genome_version} absolute_cn/wisecondorx/${bam_type}/${binsize}kb
    """
}

