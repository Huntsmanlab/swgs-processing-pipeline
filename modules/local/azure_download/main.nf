process AZURE_DOWNLOAD {
    tag "Downloading samples from Azure"
    publishDir params.input_directory, pattern:'reads', mode:'copy'
    publishDir params.input_directory, pattern:'sample_fastq.csv', mode:'copy'
    secret "AZ_SAS_TOKEN"

    container 'docker://huntsmanlab/azure_download:latest'	
    cpus 1
    memory '8 GB'
    time '2h'

    input:
    path(script)
	path(csv_file)
	val(paired_end)

    output:
    path('reads')
    path('sample_fastq_workdir.csv')
    path('sample_fastq.csv')

    script:
    """
    mkdir -p reads
	env AZCOPY_LOG_LOCATION=./.azcopy AZCOPY_JOB_PLAN_LOCATION=./.azcopy Rscript ${script} \$AZ_SAS_TOKEN ${csv_file} ${paired_end} reads ${params.input_directory}    
    """
}
