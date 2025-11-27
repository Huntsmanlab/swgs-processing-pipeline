process TRIMMOMATIC {
    tag "trimming reads of $sample_id"

	container 'docker://staphb/trimmomatic:latest'    
	cpus 8
	memory '32 GB'
	time '10h'

    input:
    tuple val(sample_id), path(reads)

    output:
	tuple val(sample_id), path("trimmed/${sample_id}/*")

    script:
	"""
	mkdir -p "trimmed/${sample_id}"

	if [ ${params.crop_50} = true ]; then  
		if [ ${params.paired_end} = true ]; then 
			trimmomatic PE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}/${sample_id}.f.paired.fastq \
                ${sample_id}.f.unpaired.fastq \
                trimmed/${sample_id}/${sample_id}.r.paired.fastq \
                ${sample_id}.r.unpaired.fastq \
                ILLUMINACLIP:/usr/local/bin/adapters/TruSeq3-PE-2.fa:2:30:10:2:true \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5 CROP:50
		else
            trimmomatic SE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}/${sample_id}.se.fastq \
                ILLUMINACLIP:/usr/local/bin/adapters/TruSeq3-SE.fa:2:30:10:2 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5 CROP:50
		fi
	else
		if [ ${params.paired_end} = true ]; then
            trimmomatic PE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}/${sample_id}.f.paired.fastq \
                ${sample_id}.f.unpaired.fastq \
                trimmed/${sample_id}/${sample_id}.r.paired.fastq \
                ${sample_id}.r.unpaired.fastq \
                ILLUMINACLIP:/usr/local/bin/adapters/TruSeq3-PE-2.fa:2:30:10:2:true \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5
		else
            trimmomatic SE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}/${sample_id}.se.fastq \
                ILLUMINACLIP:/usr/local/bin/adapters/TruSeq3-SE.fa:2:30:10:2 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5
		fi
	fi
	"""
}
