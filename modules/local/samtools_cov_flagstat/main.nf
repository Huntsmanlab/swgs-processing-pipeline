process SAMTOOLS_COV_FLAGSTAT {
    tag "generating alignment metrics for  $sample_id"
    publishDir params.output_directory, mode:'copy'
	container 'docker://staphb/samtools:1.19'
	cpus 8
	memory '32 G'
	time '2h'    

    input:
    tuple val(sample_id), path(bam), path(bai), path(txt)

    output:
	path("alignment_metrics/${sample_id}/*")
	path("alignment_metrics/${sample_id}/${sample_id}.flagstat.txt"), emit: log1
	path("alignment_metrics/${sample_id}/${sample_id}.coverageTable.tsv"), emit: log2
    
    script:
    """
    mkdir -p alignment_metrics/${sample_id}
    samtools index -@ ${task.cpus} ${bam}
    samtools coverage ${bam} > alignment_metrics/${sample_id}/${sample_id}.coverageTable.tsv
    samtools flagstat -@ ${task.cpus} ${bam} > alignment_metrics/${sample_id}/${sample_id}.flagstat.txt
    """
}
