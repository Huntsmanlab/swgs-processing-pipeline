process SAMTOOLS_SORT_INDEX {
	tag "sorting and indexing bam of $sample_id"
	container 'docker://staphb/samtools:1.19'
	cpus 8
	memory '32 G'
	time '2h'
	
	input: 
	tuple val(sample_id), path(sam)

	output:
	tuple val(sample_id), path("samtools/${sample_id}/${sample_id}.bwa.sorted.bam")

	script:
	"""
	mkdir -p "samtools/${sample_id}"

	# sam to bam piped to sort bam
	samtools view -hbS -@ ${task.cpus} ${sam} | \
		samtools sort -m 2G -@ ${task.cpus} -o samtools/${sample_id}/${sample_id}.bwa.sorted.bam

	# index bam
	samtools index -@ ${task.cpus} samtools/${sample_id}/${sample_id}.bwa.sorted.bam
	"""
}
