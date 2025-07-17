process FASTQC_PREALIGN {
    tag "pre-alignment fastqc on $sample_id"
    publishDir params.output_directory, mode:'copy', pattern: 'reports/pre_alignment_fastqc_reports/*'
    container = 'docker://staphb/fastqc:latest'
    cpus 1
    memory '8 GB'
    time '4h'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs", emit: logs
    path "reports/pre_alignment_fastqc_reports/*" 

    script:
    """
    TEMP=\$(mktemp -d --tmpdir=.)
    mkdir fastqc_${sample_id}_logs
    mkdir -p reports/pre_alignment_fastqc_reports
    env _JAVA_OPTIONS='-XX:-UsePerfData' fastqc --dir \$TEMP -o fastqc_${sample_id}_logs -f fastq ${reads.join(' ')}
    cp fastqc_${sample_id}_logs/${sample_id}_*fastqc.html reports/pre_alignment_fastqc_reports/
    """
}
