process ALIGNMENT {
    tag "Aligning ${sample_id}"
	container 'docker://clinicalgenomics/bwa-mem2:2.2.1'
	cpus 16
	memory '64 GB'
    time '18h'    

    input:
    tuple val(sample_id), path(reads)
	path(ref_genome_directory)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.bwa.sam")

    script:
    """
	mkdir -p "alignment_${sample_id}"
    bwa-mem2 mem -M -t ${task.cpus} ${params.reference_genome_path} ${reads} -R '@RG\\tID:${sample_id}_ID\\tSM:${sample_id}\\tLB:${sample_id}_LB\\tPL:ILLUMINA' -o alignment_${sample_id}/${sample_id}.bwa.sam 
    """      
}
