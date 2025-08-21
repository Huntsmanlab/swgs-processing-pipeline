process MULTIQC_POSTALIGN {
	tag "post-alignment multiqc"
    publishDir params.output_directory, mode:'copy'
	container 'docker://quay.io/biocontainers/multiqc:1.3--py35_2'
	cpus 1
	memory '8 GB'
	time '1h'

    input:
    path '*'
    path (multiqc_config)
 
    output:
    path "reports/final_multiqc_report.html"
 
    script:
    """
    mkdir -p reports
    multiqc .
    mv multiqc_report.html reports/final_multiqc_report.html
    """
}
