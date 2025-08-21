process MARK_DUPLICATES {
	tag "marking duplicates reads in $sample_id"
	publishDir params.output_directory, mode:'copy'
	container 'docker://broadinstitute/picard:latest'
	cpus 1
	memory '128 G'
	time '12h'

	input:
	tuple val(sample_id), path(bam)

	output:
	tuple val(sample_id), path("alignment_outputs/${sample_id}/${sample_id}.*bwa.sorted.mkdup.bam"), path("alignment_outputs/${sample_id}/${sample_id}.*bwa.sorted.mkdup.bai"), path("alignment_outputs/${sample_id}/${sample_id}.markdups_metrics.txt")
	path("alignment_outputs/${sample_id}/${sample_id}.*bwa.sorted.mkdup.bam")
	path("alignment_outputs/${sample_id}/${sample_id}.markdups_metrics.txt"), emit: log

	script:
	"""
	TEMP=\$(mktemp -d --tmpdir=.)
    mkdir -p "alignment_outputs/${sample_id}"

	if [ ${params.paired_end} = true ]; then

    	if [ ${params.crop_50} = true ]; then
        	java -Xmx${task.memory.toGiga() - 1}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            	--INPUT ${bam} \
            	--OUTPUT alignment_outputs/${sample_id}/${sample_id}.pe.50bp.bwa.sorted.mkdup.bam \
            	--METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.markdups_metrics.txt \
            	--CREATE_INDEX true \
				--TMP_DIR \$TEMP
    	else
        	java -Xmx${task.memory.toGiga() - 1}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
           		--INPUT ${bam} \
            	--OUTPUT alignment_outputs/${sample_id}/${sample_id}.pe.bwa.sorted.mkdup.bam \
            	--METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.markdups_metrics.txt \
            	--CREATE_INDEX true \
				--TMP_DIR \$TEMP
    	fi
	else
		if [ ${params.crop_50} = true ]; then
            java -Xmx${task.memory.toGiga() - 1}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
                --INPUT ${bam} \
                --OUTPUT alignment_outputs/${sample_id}/${sample_id}.se.50bp.bwa.sorted.mkdup.bam \
                --METRICS_FILE alignment_outputs/${sample_id}/${sample_id}markdups_metrics.txt \
                --CREATE_INDEX true \
				--TMP_DIR \$TEMP
        else
            java -Xmx${task.memory.toGiga() - 1}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
                --INPUT ${bam} \
                --OUTPUT alignment_outputs/${sample_id}/${sample_id}.se.bwa.sorted.mkdup.bam \
                --METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.markdups_metrics.txt \
                --CREATE_INDEX true \
				--TMP_DIR \$TEMP
        fi
	fi
 	"""
}
