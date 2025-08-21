process FASTQC_POSTALIGN {
    tag "post-alignment fastqc on $sample_id"
    publishDir params.output_directory, mode:'copy', pattern: 'reports/post_alignment_fastqc_reports/**'
    container = 'docker://staphb/fastqc:latest'
    cpus 1
    memory '8 GB'
    time '4h'

	input:
	tuple val(sample_id), path(bam), path(bai), path(txt)

	output:
    path "reports/post_alignment_fastqc_reports/**"
	path "post_alignment_fastqc_logs_${sample_id}", emit: log

	script:
    """
    TEMP=\$(mktemp -d --tmpdir=.)
    mkdir -p post_alignment_fastqc_logs_${sample_id}
    mkdir -p reports/post_alignment_fastqc_reports
    env _JAVA_OPTIONS='-XX:-UsePerfData' fastqc ${bam} -o post_alignment_fastqc_logs_${sample_id} --dir \$TEMP
    cp post_alignment_fastqc_logs_${sample_id}/${sample_id}*.html reports/post_alignment_fastqc_reports
    """
}
