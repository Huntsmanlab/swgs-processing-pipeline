process ICHOR_PROCESS {
    publishDir params.output_directory, mode:'copy'
    tag "Processing ichor output ${sample_id}"
    container 'docker://huntsmanlab/utanos:latest' 
    cpus 1
    memory '8 G'
    time '2h'

	input:
    tuple path(bins_file), path(segs_file), val(sample_id), val(bamtype), val(binsize)

    output:
    tuple val(bamtype), val(binsize), path("relative_cn/ichor/${bamtype}/${binsize}kb/results/${sample_id}_rCN.tsv")
    
    script:
    """
    mkdir -p "relative_cn/ichor/${bamtype}/${binsize}kb/results"
    process_ichor.R ${sample_id} ${bins_file} ${segs_file} "relative_cn/ichor/${bamtype}/${binsize}kb/results/" ${params.reference_genome_version}     
    """

}
