process AZURE_DOWNLOAD {
    tag "Downloading samples from Azure"
    publishDir params.input_directory, mode:'copy'
    secret "AZ_SAS_TOKEN"

    container 'docker://huntsmanlab/azure_download:latest'	
    cpus 1
    memory '8 GB'
    time '2h'

    input:
	path(csv_file)
	val(paired_end)

    output:
    path('reads')
    path('sample_fastq.csv')

    script:
    """
    mkdir -p reads/azure
	env AZCOPY_LOG_LOCATION=./.azcopy AZCOPY_JOB_PLAN_LOCATION=./.azcopy Rscript ${projectDir}/bin/downloadAZ.R \$AZ_SAS_TOKEN ${csv_file} ${paired_end} reads ${params.input_directory}   
    """
}
