process MARK_DUPLICATES {
	tag "marking duplicates reads in $sample_id"
	publishDir params.output_directory, mode:'copy'
	container 'docker://broadinstitute/picard:latest'
	cpus 1
	memory '32 G'
	time '12h'

	input:
	tuple val(sample_id), path(bam)

	output:
	tuple val(sample_id), 
	path("alignment_outputs/${sample_id}/${sample_id}.*bwa.sorted.mkdup.bam"),
    path("alignment_outputs/${sample_id}/${sample_id}.*bwa.sorted.mkdup.bam.bai"), 
    path("alignment_outputs/${sample_id}/${sample_id}.*mark_duplicates.txt")

	script:
	"""
    mkdir -p "alignment_outputs/${sample_id}"

	if [ ${params.paired_end} = true ]; then

    	if [ ${params.crop_50} = true ]; then
        	java -Xmx${task.memory.toGiga()}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            	--INPUT ${bam} \
            	--OUTPUT alignment_outputs/${sample_id}/${sample_id}.pe.50bp.bwa.sorted.mkdup.bam \
            	--METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.50bp.mark_duplicates.txt \
            	--CREATE_INDEX true
			mv alignment_outputs/${sample_id}/${sample_id}.pe.50bp.bwa.sorted.mkdup.bai alignment_outputs/${sample_id}/${sample_id}.pe.50bp.bwa.sorted.mkdup.bam.bai
    	else
        	java -Xmx${task.memory.toGiga()}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
           		--INPUT ${bam} \
            	--OUTPUT alignment_outputs/${sample_id}/${sample_id}.pe.bwa.sorted.mkdup.bam \
            	--METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.mark_duplicates.txt \
            	--CREATE_INDEX true
			mv alignment_outputs/${sample_id}/${sample_id}.pe.bwa.sorted.mkdup.bai alignment_outputs/${sample_id}/${sample_id}.pe.bwa.sorted.mkdup.bam.bai
    	fi
	else
		if [ ${params.crop_50} = true ]; then
            java -Xmx${task.memory.toGiga()}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
                --INPUT ${bam} \
                --OUTPUT alignment_outputs/${sample_id}/${sample_id}.se.50bp.bwa.sorted.mkdup.bam \
                --METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.50bp.mark_duplicates.txt \
                --CREATE_INDEX true
            mv alignment_outputs/${sample_id}/${sample_id}.se.50bp.bwa.sorted.mkdup.bai alignment_outputs/${sample_id}/${sample_id}.se.50bp.bwa.sorted.mkdup.bam.bai
        else
            java -Xmx${task.memory.toGiga()}g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
                --INPUT ${bam} \
                --OUTPUT alignment_outputs/${sample_id}/${sample_id}.se.bwa.sorted.mkdup.bam \
                --METRICS_FILE alignment_outputs/${sample_id}/${sample_id}.mark_duplicates.txt \
                --CREATE_INDEX true
            mv alignment_outputs/${sample_id}/${sample_id}.se.bwa.sorted.mkdup.bai alignment_outputs/${sample_id}/${sample_id}.se.bwa.sorted.mkdup.bam.bai
        fi
	fi
 	"""
}
