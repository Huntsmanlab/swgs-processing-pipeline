process MULTIQC_PREALIGN {
    publishDir params.output_directory, mode:'copy'
	tag "pre-alignment multiqc"	

	container 'docker://quay.io/biocontainers/multiqc:1.3--py35_2'
	cpus 1
	memory '8 GB'
	time '1h'

    input:
    path '*'
 
    output:
    path "reports/pre_alignment_multiqc_report.html"
 
    script:
    """
    mkdir -p reports
    multiqc .
    mv multiqc_report.html reports/pre_alignment_multiqc_report.html
    """
}
