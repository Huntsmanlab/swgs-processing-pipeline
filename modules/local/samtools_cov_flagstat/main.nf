process SAMTOOLS_COV_FLAGSTAT {
    tag "generating alignment metrics for  $sample_id"
    publishDir params.output_directory, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(metrics)

    output:
    tuple val(sample_id), 
	path("alignment_metrics/${sample_id}/${sample_id}.coverageTable.tsv"), 
    path("alignment_metrics/${sample_id}/${sample_id}.flagstat.txt")
    
    script:
    """
    mkdir -p alignment_metrics/${sample_id}

    samtools index -@ ${task.cpus} ${bam}
    samtools coverage ${bam} > alignment_metrics/${sample_id}/${sample_id}.coverageTable.tsv
    samtools flagstat -@ ${task.cpus} ${bam} > alignment_metrics/${sample_id}/${sample_id}.flagstat.txt
    """
}
