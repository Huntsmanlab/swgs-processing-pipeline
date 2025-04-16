#!/usr/bin/env nextflow

/*
 * Copyright (c) 2023, Maxwell Douglas
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

/* 
 * 'Shallow WGS - NF' - A Nextflow pipeline for processing shallow WGS data
 * 
 * This pipeline reproduces the analysis pipeline described in " Targeted 
 * and shallow whole genome sequencing identifies therapeutic opportunities 
 * in p53abn endometrial cancers"
 * 
 * Maxwell Douglas 
 */


log.info """\

S H A L L O W   W G S   C N - C A L L I N G  -  N F    v 2.1 
================================
genome   : $params.ref_path
reads    : $params.reads
results  : $params.outdir
"""


// Download fastqs from Azure
process AZURE_DOWNLOAD {
    tag "Download samples from Azure"
    publishDir params.indir, pattern:'reads', mode:'copy'
    publishDir params.indir, pattern:'sample_fastq_pub.csv', mode:'copy'

    secret "AZ_SAS_TOKEN"

    input:
    path(script)

    output:
    path('reads')
    path('sample_fastq.csv')
    path('sample_fastq_pub.csv')

    script:
    """
    mkdir -p reads
    env AZCOPY_LOG_LOCATION=./.azcopy AZCOPY_JOB_PLAN_LOCATION=./.azcopy Rscript ${script} \$AZ_SAS_TOKEN ${params.az_csv} ${params.pairedend} reads ${params.indir}
    """
}

// Get reference genome
process DOWNLOAD_MM10 {
    errorStrategy 'retry', maxRetries: 2
    tag "Download Reference Genome mm10"
    
    output:
    path "reference_genome/mm10.fa.gz"

    script:
    """
    mkdir -p reference_genome
    curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz > reference_genome/mm10.fa.gz
    """
}
process DOWNLOAD_HG19 {
    errorStrategy 'retry', maxRetries: 2
    tag "Download Reference Genome hg19"
    
    output:
    path "reference_genome/hg19.fa.gz"

    script:
    """
    mkdir -p reference_genome
    curl https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz > reference_genome/hg19.fa.gz
    """
}
process DOWNLOAD_HG38 {
    errorStrategy 'retry', maxRetries: 2
    tag "Download Reference Genome hg38"
    
    output:
    path "reference_genome/hg38.fa.gz"

    script:
    """
    mkdir -p reference_genome
    curl https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz > reference_genome/hg38.fa.gz
    """
}

// Pre-alignment QC
process FASTQC1 {
    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode:'copy', pattern: 'reports/pre_alignment_fastqc_reports/*'

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
    env _JAVA_OPTIONS='-XX:-UsePerfData' fastqc --dir \$TEMP -o fastqc_${sample_id}_logs -f fastq $reads
    cp fastqc_${sample_id}_logs/${sample_id}_*fastqc.html reports/pre_alignment_fastqc_reports/
    """
}

// Create a QC report
process MULTIQC1 {
    publishDir params.outdir, mode:'copy'

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

// Index reference genome for bwa-mem2
process INDEX_BWAREF {
    tag "Index Reference Genome $genome"
    publishDir params.outdir, mode:'copy'

    input:
    file genome

    output:
    path "${genome.simpleName}_bwa_refgenome_index/${genome.simpleName}.*"

    script:
    """
    mkdir -p ${genome.simpleName}_bwa_refgenome_index
    cp ${genome} ${genome.simpleName}_bwa_refgenome_index/${genome}
    bwa-mem2 index ${genome} -p ${genome.simpleName}_bwa_refgenome_index/${genome}
    """
}

// Trim reads using Trimmomatic (paired-end)
process TRIMMOMATIC_PE {
    tag "TRIMMOMATIC_PE on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed/${sample_id}.{f,r}.paired.fastq"), 
    path("trimmed/${sample_id}.{f,r}.unpaired.fastq")

    script:
        """
        mkdir trimmed

        if [ ${params.crop50} = true ]; then
            trimmomatic PE -threads ${task.ncpus} -phred33 $reads \
                trimmed/${sample_id}.f.paired.fastq \
                trimmed/${sample_id}.f.unpaired.fastq \
                trimmed/${sample_id}.r.paired.fastq \
                trimmed/${sample_id}.r.unpaired.fastq \
                ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:true \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5 CROP:50
        else
            trimmomatic PE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}.f.paired.fastq \
                trimmed/${sample_id}.f.unpaired.fastq \
                trimmed/${sample_id}.r.paired.fastq \
                trimmed/${sample_id}.r.unpaired.fastq \
                ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:true \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5
        fi
        """
}

// Trim reads using Trimmomatic (single-end)
process TRIMMOMATIC_SE {
    tag "TRIMMOMATIC_SE on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed/${sample_id}.se.fastq")

    script:
        """
        mkdir trimmed

        if [ ${params.crop50} = true ]; then
            trimmomatic SE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}.se.fastq \
                ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:2 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5 CROP:50
        else 
            trimmomatic SE -threads ${task.cpus} -phred33 $reads \
                trimmed/${sample_id}.se.fastq \
                ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:2 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5
        fi
        """
}

// minimap2 alignment process
process MINIMAP2_ALIGNMENT {
    tag "MINIMAP2 ALIGNMENT on $sample_id"
    
    input:
    tuple val(sample_id), path(reads), path(unpaired_reads)
    path(reference)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam")

    script:
    """
    mkdir -p "alignment_${sample_id}"
    minimap2 -a -x sr -Y -K 100M -t ${task.cpus} ${reference} ${reads} | samtools view -hbS | \
        samtools sort -m 2G -@ ${task.cpus} -o alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam
    samtools index -@ ${task.cpus} alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam
    """
}

// bwa alignment process (pe)
process BWAMEM2_ALIGNMENT_PE {
    tag "BWA-MEM 2 ALIGNMENT (pe) on $sample_id"
    
    input:
    tuple val(sample_id), path(paired_reads), path(unpaired_reads)
    path(bwa_reference)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.pe.bwa.sam"),
    path("alignment_${sample_id}/${sample_id}.up.bwa.sam")

    script:
    // bwaList = bwa_reference.collect()
    // fa_path = bwaList[4]                           // Access individual elements from the list
    unpairedList = unpaired_reads.collect()
    forward_unpaired = unpairedList[0]
    """
    mkdir -p "alignment_${sample_id}"
    bwa-mem2 mem -M -t ${task.cpus} ${params.ref_path} ${paired_reads} -R '@RG\\tID:${sample_id}_ID\\tSM:${sample_id}\\tLB:${sample_id}_LB\\tPL:ILLUMINA' > alignment_${sample_id}/${sample_id}.pe.bwa.sam 

    bwa-mem2 mem -M -t ${task.cpus} ${params.ref_path} ${forward_unpaired} -R '@RG\\tID:${sample_id}_ID\\tSM:${sample_id}\\tLB:${sample_id}_LB\\tPL:ILLUMINA' > alignment_${sample_id}/${sample_id}.up.bwa.sam
    """      
}//.after(INDEX_BWAREF)

// bwa alignment process (se)
process BWAMEM2_ALIGNMENT_SE {
    tag "BWA-MEM 2 ALIGNMENT (se) on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)
    path(bwa_reference)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.se.bwa.sam")

    script:
    // bwaList = bwa_reference.collect()
    // fa_path = bwaList[4]                           // Access individual elements from the list
    """
    mkdir -p "alignment_${sample_id}"
    bwa-mem2 mem -M -t ${task.cpus} ${params.ref_path} ${reads} -R '@RG\\tID:${sample_id}_ID\\tSM:${sample_id}\\tLB:${sample_id}_LB\\tPL:ILLUMINA' > alignment_${sample_id}/${sample_id}.se.bwa.sam
    """
}//.after(INDEX_BWAREF)

// sorting and indexing bam (pe)
process SORT_INDEX_PE {
    tag "SORTING AND INDEXING (pe) on $sample_id"
    
    input:
    tuple val(sample_id), path(sam_pe), path(sam_up) 

    output:
    tuple val(sample_id), path("sorting_${sample_id}/${sample_id}.pe.bwa.sorted.bam"), path("sorting_${sample_id}/${sample_id}.se.merged.sorted.bam")

    script:
    """
    mkdir -p "sorting_${sample_id}"

    samtools view -hbS ${sample_id}.pe.bwa.sam | \
        samtools sort -m 2G -@ ${task.cpus} -o sorting_${sample_id}/${sample_id}.pe.bwa.sorted.bam
    samtools index -@ ${task.cpus} sorting_${sample_id}/${sample_id}.pe.bwa.sorted.bam

    samtools view -hbS ${sample_id}.up.bwa.sam | \
        samtools sort -m 2G -@ ${task.cpus} -o sorting_${sample_id}/${sample_id}.up.bwa.sorted.bam
    samtools index -@ ${task.cpus} sorting_${sample_id}/${sample_id}.up.bwa.sorted.bam

    # Grab forward strand from P.E. post_alignment
    samtools view -F 16 -o sorting_${sample_id}/${sample_id}.pe.f.bwa.sorted.bam \
        sorting_${sample_id}/${sample_id}.pe.bwa.sorted.bam
    samtools merge -f -@ ${task.cpus} sorting_${sample_id}/${sample_id}.se.merged.bam \
        sorting_${sample_id}/${sample_id}.up.bwa.sorted.bam \
        sorting_${sample_id}/${sample_id}.pe.f.bwa.sorted.bam
    # Re-Sort now merged BAM
    samtools sort -m 2G -@ ${task.cpus} -o sorting_${sample_id}/${sample_id}.se.merged.sorted.bam \
        sorting_${sample_id}/${sample_id}.se.merged.bam
    samtools index -@ ${task.cpus} sorting_${sample_id}/${sample_id}.se.merged.sorted.bam
    """
}//.after(INDEX_BWAREF)

// sorting and indexing bam (se)
process SORT_INDEX_SE {
    tag "SORTING AND INDEXING (se) on $sample_id"
    
    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("sorting_${sample_id}/${sample_id}.se.bwa.sorted.bam")

    script:
    """
    mkdir -p "sorting_${sample_id}"
    samtools view -hbS ${sample_id}.se.bwa.sam | \
        samtools sort -m 2G -@ ${task.cpus} -o sorting_${sample_id}/${sample_id}.se.bwa.sorted.bam
    samtools index -@ ${task.cpus} sorting_${sample_id}/${sample_id}.se.bwa.sorted.bam
    """
}//.after(INDEX_BWAREF)

// Mark duplicate reads using picardtools (pe)
process MARKDUP_PE {
    tag "MARKDUP (pe) on $sample_id"
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam_pe), path(bam_se)

    output:
    tuple val(sample_id), 
    path("processing_outputs/${sample_id}/${sample_id}.pe.*bwa.sorted.mkdup.bam"), 
    path("processing_outputs/${sample_id}/${sample_id}.pe.*bwa.sorted.mkdup.bai"), 
    path("processing_outputs/${sample_id}/${sample_id}.se.*merged.sorted.mkdup.bam"), 
    path("processing_outputs/${sample_id}/${sample_id}.se.*merged.sorted.mkdup.bai"), 
    path("processing_outputs/${sample_id}/${sample_id}.pe.*marked_duplicates.metrics.txt"),
    path("processing_outputs/${sample_id}/${sample_id}.se.*marked_duplicates.metrics.txt")
    
    script:
    """
    TEMP=\$(mktemp -d --tmpdir=.)
    mkdir -p "processing_outputs/${sample_id}"

    if [ ${params.crop50} = true ]; then
        java -Xmx16g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            -I ${bam_pe} \
            -O processing_outputs/${sample_id}/${sample_id}.pe.50bp.bwa.sorted.mkdup.bam \
            -M processing_outputs/${sample_id}/${sample_id}.pe.50bp.marked_duplicates.metrics.txt \
            --CREATE_INDEX true \
            --TMP_DIR \$TEMP
        java -Xmx16g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            -I ${bam_se} \
            -O processing_outputs/${sample_id}/${sample_id}.se.50bp.merged.sorted.mkdup.bam \
            -M processing_outputs/${sample_id}/${sample_id}.se.50bp.marked_duplicates.metrics.txt \
            --CREATE_INDEX true \
            --TMP_DIR \$TEMP
    else
        java -Xmx16g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            -I ${bam_pe} \
            -O processing_outputs/${sample_id}/${sample_id}.pe.bwa.sorted.mkdup.bam \
            -M processing_outputs/${sample_id}/${sample_id}.pe.marked_duplicates.metrics.txt \
            --CREATE_INDEX true \
            --TMP_DIR \$TEMP
        java -Xmx16g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            -I ${bam_se} \
            -O processing_outputs/${sample_id}/${sample_id}.se.merged.sorted.mkdup.bam \
            -M processing_outputs/${sample_id}/${sample_id}.se.marked_duplicates.metrics.txt \
            --CREATE_INDEX true \
            --TMP_DIR \$TEMP
    fi
    """
}

// Mark duplicate reads using picardtools (se)
process MARKDUP_SE {
    tag "MARKDUP (se) on $sample_id"
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("processing_outputs/${sample_id}/${sample_id}.se.*bwa.sorted.mkdup.bam"), 
    path("processing_outputs/${sample_id}/${sample_id}.se.*bwa.sorted.mkdup.bai"), 
    path("processing_outputs/${sample_id}/${sample_id}.se.*marked_duplicates.metrics.txt")
    
    script:
    """
    mkdir -p "processing_outputs/${sample_id}"

    if [ ${params.crop50} = true ]; then
        java -Xmx16g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            I=${bam} \
            O=processing_outputs/${sample_id}/${sample_id}.se.50bp.bwa.sorted.mkdup.bam \
            M=processing_outputs/${sample_id}/${sample_id}.se.50bp.marked_duplicates.metrics.txt \
            CREATE_INDEX=true
    else
        java -Xmx16g -XX:-UsePerfData -jar /usr/picard/picard.jar MarkDuplicates \
            I=${bam} \
            O=processing_outputs/${sample_id}/${sample_id}.se.bwa.sorted.mkdup.bam \
            M=processing_outputs/${sample_id}/${sample_id}.se.marked_duplicates.metrics.txt \
            CREATE_INDEX=true
    fi
    """
}

// Calculating a few samtools alignment metrics (pe)
process INDCOVFLAG_PE {
    tag "Re-index and sort (pe) on $sample_id"
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam_pe), path(bai_pe), path(bam_se), path(bai_se), path(metrics_pe), path(metrics_se)

    output:
    tuple val(sample_id), 
    path("metrics_files/${sample_id}.pe.coverageTable.tsv"), 
    path("metrics_files/${sample_id}.pe.flagstat.txt"), 
    path("metrics_files/${sample_id}.se.coverageTable.tsv"), 
    path("metrics_files/${sample_id}.se.flagstat.txt")
    
    script:
    """
    mkdir -p metrics_files

    samtools index -@ ${task.cpus} ${bam_pe}
    samtools coverage ${bam_pe} > metrics_files/${sample_id}.pe.coverageTable.tsv
    samtools flagstat -@ ${task.cpus} ${bam_pe} > metrics_files/${sample_id}.pe.flagstat.txt

    samtools index -@ ${task.cpus} ${bam_se}
    samtools coverage ${bam_se} > metrics_files/${sample_id}.se.coverageTable.tsv
    samtools flagstat -@ ${task.cpus} ${bam_se} > metrics_files/${sample_id}.se.flagstat.txt
    """
}

// Calculating a few samtools alignment metrics (se)
process INDCOVFLAG_SE {
    tag "Re-index and sort (se) on $sample_id"
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(metrics)

    output:
    tuple val(sample_id), path("metrics_files/${sample_id}.se.coverageTable.tsv"), 
    path("metrics_files/${sample_id}.se.flagstat.txt")
    
    script:
    """
    mkdir -p metrics_files
    samtools index -@ ${task.cpus} ${bam}
    samtools coverage ${bam} > metrics_files/${sample_id}.se.coverageTable.tsv
    samtools flagstat -@ ${task.cpus} ${bam} > metrics_files/${sample_id}.se.flagstat.txt
    """
}

// Post-alignment QC (pe)
process FASTQC2_PE {
    tag "FASTQC (pe) on $sample_id"
    publishDir params.outdir, mode: 'copy', pattern: 'reports/post_alignment_fastqc_reports/*'

    input:
    tuple val(sample_id), path(bam_pe), path(bai_pe), path(bam_se), path(bai_se), path(metrics_pe), path(metrics_se)
 
    output:
    path "post_alignment_fastqc_logs_${sample_id}", emit: logs
    path "reports/post_alignment_fastqc_reports/*"
 
    script:
    """
    TEMP=\$(mktemp -d --tmpdir=.)
    mkdir -p post_alignment_fastqc_logs_${sample_id}
    mkdir -p reports/post_alignment_fastqc_reports
    env _JAVA_OPTIONS='-XX:-UsePerfData' fastqc ${bam_pe} ${bam_se} -o post_alignment_fastqc_logs_${sample_id} --dir \$TEMP
    cp post_alignment_fastqc_logs_${sample_id}/${sample_id}*.html reports/post_alignment_fastqc_reports
    """
}

// Post-alignment QC (se)
process FASTQC2_SE {
    tag "FASTQC (se) on $sample_id"
    publishDir params.outdir,  mode:'copy', pattern: 'reports/post_alignment_fastqc_reports/*'

    input:
    tuple val(sample_id), path(bam), path(bai), path(metrics)
 
    output:
    path "post_alignment_fastqc_logs_${sample_id}", emit: logs
    path "reports/post_alignment_fastqc_reports/*"
 
    script:
    """
    TEMP=\$(mktemp -d --tmpdir=.)
    mkdir -p post_alignment_fastqc_logs_${sample_id}
    mkdir -p reports/post_alignment_fastqc_reports
    env _JAVA_OPTIONS='-XX:-UsePerfData' fastqc ${bam} -o post_alignment_fastqc_logs_${sample_id} --dir \$TEMP
    cp post_alignment_fastqc_logs_${sample_id}/${sample_id}*.html reports/post_alignment_fastqc_reports
    """
}

// Create a QC report
process MULTIQC2 {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'
    path (multiqc_config)
 
    output:
    path "reports/final_multiqc_report.html"
 
    script:
    """
    mkdir -p reports
    multiqc .
    mv multiqc_report.html reports/final_multiqc_report.html
    """
}

// Bin annotation for QDNAseq
process QDNA_BINS {
    input:
    path(script)
    val(binsize)
    path(bwavgbed)

    output:
    tuple val(binsize), path("qd_bins")

    script: 
    """
        mkdir -p qd_bins
        CPATH=\$(pwd)
        Rscript ${script} ${task.cpus} ${binsize} ${params.genome} ${params.qd_mappability} ${params.qd_blacklist} \$CPATH/${bwavgbed} ${params.qd_nbams} qd_bins ${params.pairedend}
    """
}

// Call Copy-Numbers using QDNAseq
process CN_QDNA1 {
    publishDir params.outdir, mode:'copy'

    input:
    path(script)
    tuple val(binsize), path(binannos_dir), val(bam_type), val(bams)

    output:
    path("relative_cns/qdnaseq/${bam_type}/${binsize}kb")
    val(binsize)
    val(bam_type)

    script:
    """
    mkdir -p "relative_cns/qdnaseq/${bam_type}/${binsize}kb"
    printf '%s\n' "${bams.join('\n')}" > bamfileslist.txt
    shopt -s nocaseglob extglob
    Rscript ${script} ${binsize}kb ${task.cpus} relative_cns/qdnaseq/${bam_type}/${binsize}kb bamfileslist.txt ${binannos_dir}/*@(${binsize}kb|${bam_type})*@(${binsize}kb|${bam_type})*.rds
    rm bamfileslist.txt
    """
}

// Reference creation for WisecondorX
process WX_REF {
    input:
    path(normals)
    val(binsize)
    path(bams)

    output:
    tuple val(binsize), path("wx_references")

    script: 
    """
        mkdir -p wx_references

        if [ ${params.wx_newref_frombam} = true ]; then
            mkdir -p converts 
            for FILE in ${bams}/*.bam; do

                SAMPLE=\$(basename "\${FILE%%.*}")
                if [[ \$FILE == *".pe."* ]]; then
                    WisecondorX convert \$FILE converts/\${SAMPLE}.pe.npz --binsize ${binsize}000
                else
                    WisecondorX convert \$FILE converts/\${SAMPLE}.se.npz --binsize ${binsize}000
                fi

            done
            NORMALS="converts"
        else
            NORMALS="${normals}"
        fi

        WisecondorX newref \$NORMALS/*.se.npz wx_references/reference_${binsize}kb.se.npz \
        --binsize ${binsize}000 --cpus ${task.cpus} --yfrac 1

        if [ ${params.pairedend} = true ]; then
            WisecondorX newref \$NORMALS/*.pe.npz wx_references/reference_${binsize}kb.pe.npz \
            --binsize ${binsize}000 --cpus ${task.cpus} --yfrac 1
        fi
    """
}

// Call Copy-Numbers using WisecondorX
process CN_WX1 {
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(binsize), path(wx_ref), val(bam_type), val(bams)

    output:
    path("relative_cns/wisecondorx/${bam_type}/${binsize}kb")
    val(binsize)
    val(bam_type)

    script:
    """
    mkdir -p "relative_cns/wisecondorx/${bam_type}/${binsize}kb"
    mkdir -p converts
    printf '%s\n' "${bams.join('\n')}" > bamfileslist.txt

    while read LINE; do

        SAMPLE=\$(basename "\${LINE%%.*}")

        WisecondorX convert \$LINE converts/\${SAMPLE}.${bam_type}.npz --binsize ${binsize}000

    done < bamfileslist.txt

    for FILE in converts/*; do

        SAMPLE=\$(basename "\${FILE%%.*}")

        mkdir -p relative_cns/wisecondorx/${bam_type}/${binsize}kb/\${SAMPLE}

        if [[ \$FILE == *"pe."* ]]; then
            WisecondorX predict \$FILE ${wx_ref}/*${binsize}kb*.pe.npz relative_cns/wisecondorx/${bam_type}/${binsize}kb/\${SAMPLE}/\${SAMPLE} --gender F --plot --bed
        else
            WisecondorX predict \$FILE ${wx_ref}/*${binsize}kb*.se.npz relative_cns/wisecondorx/${bam_type}/${binsize}kb/\${SAMPLE}/\${SAMPLE} --gender F --plot --bed
        fi

    done

    rm bamfileslist.txt
    """
}

// Scale relative copy-number to absolute copy-number from QDNAseq output
process RCN_TO_ACN_QDNA {
    publishDir params.outdir, mode:'copy'

    input:
    path(script)
    path(rcn_dir)
    val(binsize)
    val(bam_type)

    output:
    path("absolute_cns/qdnaseq/${bam_type}/${binsize}kb")

    script:
    """
    mkdir -p "absolute_cns/qdnaseq/${bam_type}/${binsize}kb"
    Rscript ${script} ${task.cpus} qdnaseq ${rcn_dir} ${params.genome} absolute_cns/qdnaseq/${bam_type}/${binsize}kb
    """
}

// Scale relative copy-number to absolute copy-number from WisecondorX output
process RCN_TO_ACN_WX {
    publishDir params.outdir, mode:'copy'

    input:
    path(script)
    path(rcn_dir)
    val(binsize)
    val(bam_type)

    output:
    path("absolute_cns/wisecondorx/${bam_type}/${binsize}kb")

    script:
    """
    mkdir wx_files
    cp ${rcn_dir}/*/*.txt wx_files && cp ${rcn_dir}/*/*.bed wx_files
    mkdir -p "absolute_cns/wisecondorx/${bam_type}/${binsize}kb"
    Rscript ${script} ${task.cpus} wx wx_files ${params.genome} absolute_cns/wisecondorx/${bam_type}/${binsize}kb
    """
}

process ICHOR_CNA {
    tag "ichorCNA ($bamtype) on $sample_id with bin size $binsize kb"
    publishDir params.outdir, mode:'copy'   
 
    input: 
    tuple val(binsize), val(sample_id), val(bamtype), path(bam), path(bamindex)

    output:
    tuple path("relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.correctedDepth.txt"), path("relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.seg"), val(sample_id), val(bamtype), val(binsize), emit: results
    path "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.cna.seg"
    path "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.correctedDepth.txt"
    path "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.seg"
    path "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.seg.txt"
    path "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/**/*"
 
    script:
    def binsize_bases = binsize.toInteger() * 1000
    def centromere_file
    if ("${params.genome}" == "hg19" || "${params.genome}" == "GRCh37") {
	centromere_file = "GRCh37.p13_centromere_UCSC-gapTable.txt"
    } else {
	centromere_file = "GRCh38.GCA_000001405.2_centromere_acen.txt"
    } 
    """
    mkdir -p "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/"
    bamindex_file=`readlink -f "${bamindex}"`
    bam_file=`readlink -f "${bam}"`
    flock --verbose --wait 30 -E 3 "\${bamindex_file}.lock" cp -n \${bamindex_file} \${bamindex_file%.bai}.bam.bai    

    readCounter --window "${binsize_bases}" --quality 20 \
    --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" \
    \${bam_file} > "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.readcounts.wig"
    
    Rscript /usr/local/bin/runIchorCNA.R \
    --id "${sample_id}" \
    --WIG "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}/${sample_id}.readcounts.wig" \
    --gcWig "/usr/local/bin/ichorCNA/inst/extdata/gc_${params.genome}_${binsize}kb.wig" \
    --mapWig "/usr/local/bin/ichorCNA/inst/extdata/map_${params.genome}_${binsize}kb.wig" \
    --centromere "/usr/local/bin/ichorCNA/inst/extdata/${centromere_file}" \
    --ploidy "${params.ichor_ploidy}"  --normal "${params.ichor_normal}" \
    --maxCN "${params.ichor_maxCN}" \
    --minMapScore 0.75 \
    --includeHOMD True \
    --estimateScPrevalence False --txnE 0.9999 --txnStrength 10000 \
    --fracReadsInChrYForMale 0.002 \
    --outDir "relative_cns/ichorCNA/${bamtype}/${binsize}kb/${sample_id}"
    """
}

process PROCESS_ICHOR {
    tag "Process ichor ($bamtype) on $sample_id"
    publishDir params.outdir, mode:'copy'

    input: 
    path(script)
    tuple path(bins_file), path(segs_file), val(sample_id), val(bamtype), val(binsize)

    output:
    tuple val(bamtype), val(binsize), path("relative_cns/ichorCNA/${bamtype}/${binsize}kb/results/${sample_id}_rCN.tsv")
    
    script:
    """
    mkdir -p "relative_cns/ichorCNA/${bamtype}/${binsize}kb/results"
    Rscript ${script} ${sample_id} ${bins_file} ${segs_file} "relative_cns/ichorCNA/${bamtype}/${binsize}kb/results/" ${params.genome}     
    """
}

process COMBINE_ICHOR {
    tag "Combine ichor $bamtype $binsize"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(bamtype), val(binsize), path(tsv_files)

    output:
    path "relative_cns/ichorCNA/combined_results/ichor_combined_${binsize}kb_${bamtype}_rCN.tsv"

    script:
    """
    mkdir -p "relative_cns/ichorCNA/combined_results/"
    head -n 1 ${tsv_files[0]} > "relative_cns/ichorCNA/combined_results/ichor_combined_${binsize}kb_${bamtype}_rCN.tsv"
    for FILE in ${tsv_files}; do
	tail -n +2 \$FILE >> "relative_cns/ichorCNA/combined_results/combined_${binsize}kb_${bamtype}_rCN.tsv"
    done
    """
}

workflow {
    // Make some channels for the needed params, reads files, and scripts
    binsizes_ch = Channel.from(params.binsizes)
    qdnaseq_script_ch = file("$projectDir/scripts/runQDNAseq.R")
    acn_script_ch = file("$projectDir/scripts/runACN.R")
    qdna_bins_script_ch = file("$projectDir/scripts/gen_bin_annot.R")
    azure_script_ch = file("$projectDir/scripts/downloadAZ.R")
    ichor_script_ch = file("$projectDir/scripts/process_ichor.R")

    if (params.from_azure) {
        az_ch = AZURE_DOWNLOAD(azure_script_ch)
        reads_ch = az_ch[1].splitCsv(skip: 1).map { row -> [row[0], row[1..-1]] }
    } else {
        if (params.use_csv) {
            match_ch = Channel.fromPath(params.samples_csv, checkIfExists: true)
            reads_ch = match_ch.splitCsv(skip: 1).map { row -> [row[0], row[1..-1]] }
        } else {
            reads_ch = Channel
                .fromFilePairs(params.reads, checkIfExists: true, size:-1) { 
                    file -> file.name.replaceAll(params.rm_regex, "")
                }
        }
    }
    
    // Get the reference genome
    // if (params.genome == 'hg19') {
    //     reference_genome_ch = DOWNLOAD_HG19()
    // } else if ( params.genome == 'hg38') {
    //     reference_genome_ch = DOWNLOAD_HG38()
    // } else if (params.genome == 'mm10') {
    //     reference_genome_ch = DOWNLOAD_MM10()
    // } else {
    //     fail "Reference genome is not recognized, please correct input string."
    // }
    
    // reads_ch.view()
    FASTQC1(reads_ch)
    // fastqc1_ch.view()
    MULTIQC1(FASTQC1.out.logs.collect())

    if (params.pairedend) {
        trimmed_ch = TRIMMOMATIC_PE(reads_ch)
    } else {
        trimmed_ch = TRIMMOMATIC_SE(reads_ch)
    }

    if (params.aligner == 'bwamem2') {
        // bwa_reference_genome_ch = INDEX_BWAREF(reference_genome_ch)
        // bwa_reference_genome_ch.view()
        if (params.pairedend) {
            align_ch = BWAMEM2_ALIGNMENT_PE(trimmed_ch, params.ref_path)
        } else {
            align_ch = BWAMEM2_ALIGNMENT_SE(trimmed_ch, params.ref_path)
        }
    } else if (params.aligner == 'minimap2') {
        align_ch = MINIMAP2_ALIGNMENT(trimmed_ch, reference_genome_ch)
    } else {
        fail "Aligner is not recognized, please correct input string. Options are 'bwamem2' or 'minimap2'."
    }

    if (params.pairedend) {
        sort_ch = SORT_INDEX_PE(align_ch)
        markdup_ch = MARKDUP_PE(sort_ch)
        indcovflag_ch = INDCOVFLAG_PE(markdup_ch)
        // fastqc2_ch = FASTQC2_PE(markdup_ch)
        FASTQC2_PE(markdup_ch)
    } else {
        sort_ch = SORT_INDEX_SE(align_ch)
        markdup_ch = MARKDUP_SE(sort_ch)
        indcovflag_ch = INDCOVFLAG_SE(markdup_ch)
        // fastqc2_ch = FASTQC2_SE(markdup_ch)
        FASTQC2_SE(markdup_ch)
    }

    if (params.runqdnaseq) {
        if (params.pairedend) {
            bamlist = markdup_ch.map { tuple -> [["pe", tuple[1]], ["se", tuple[3]]] }.flatMap().groupTuple()
	} else {
            bamlist = markdup_ch.map { tuple -> ["se", tuple[1]] }.groupTuple()
        }
        if (params.qd_new_annot) {
            qd_annot_ch = QDNA_BINS(qdna_bins_script_ch, binsizes_ch, params.qd_bwgavgbed)
            bams_ch = qd_annot_ch.combine(bamlist)
            qdna_ch = CN_QDNA1(qdnaseq_script_ch, bams_ch)
            RCN_TO_ACN_QDNA(acn_script_ch, qdna_ch)
        } else {
            bams_ch = binsizes_ch.combine(Channel.of(params.binannos)).combine(bamlist)
	    qdna_ch = CN_QDNA1(qdnaseq_script_ch, bams_ch)
	    RCN_TO_ACN_QDNA(acn_script_ch, qdna_ch)
        }
    }

    if (params.runwisex) {
        if (params.pairedend) {
            bamlist = markdup_ch.map { tuple -> [["pe", tuple[1]], ["se", tuple[3]]] }.flatMap().groupTuple()
        } else {
            bamlist = markdup_ch.map { tuple -> ["se", tuple[1]] }.groupTuple()
        }
        if (params.wx_newref) {
            wx_ref_ch = WX_REF(params.wx_normals, binsizes_ch, params.wx_nbams)
            bams_ch = wx_ref_ch.combine(bamlist)
            wx_ch = CN_WX1(bams_ch)
            RCN_TO_ACN_WX(acn_script_ch, wx_ch)
        } else {
            bams_ch = binsizes_ch.combine(Channel.of(params.wx_refs)).combine(bamlist)
            wx_ch  = CN_WX1(bams_ch)
            RCN_TO_ACN_WX(acn_script_ch, wx_ch)
        }
    }

    if (params.runichor) {
	if (params.pairedend) {
            bamlist = markdup_ch.map { vals -> [[vals[0], "pe", vals[1], vals[2]], [vals[0], "se", vals[3], vals[4]]] }.flatMap()
	    bams_ch = binsizes_ch.combine(bamlist).groupTuple(by: 3).transpose()
            bams_ch_filtered = bams_ch.filter { val -> val[0] in [10, 50, 500, 1000] } 
	    ICHOR_CNA(bams_ch_filtered)
	    PROCESS_ICHOR(ichor_script_ch, ICHOR_CNA.out.results) 
	    tsv_ch = PROCESS_ICHOR.out.groupTuple(by: [0, 1])
            COMBINE_ICHOR(tsv_ch)
	    
        } else {
            bamlist = markdup_ch.map { vals -> [vals[0], "se", vals[1], vals[2]] } 
	    bams_ch = binsizes_ch.combine(bamlist).groupTuple(by: 3).transpose()
            bams_ch_filtered = bams_ch.filter { val -> val[0] in [10, 50, 500, 1000] }
            ICHOR_CNA(bams_ch_filtered)	    
            PROCESS_ICHOR(ichor_script_ch, ICHOR_CNA.out.results)
	    tsv_ch = PROCESS_ICHOR.out.groupTuple(by: [0, 1])
	    COMBINE_ICHOR(tsv_ch)
        }
    }

    // Extract the file path from each tuple in the channel
    if (params.pairedend) {
        metrics = markdup_ch.map { tuple -> [tuple[5], tuple[6]] }
        flagstats = indcovflag_ch.map { tuple -> [tuple[2], tuple[4]] }
        combined_ch = metrics.mix(flagstats, FASTQC2_PE.out.logs).collect()
    } else {
        metrics = markdup_ch.map { tuple -> tuple[3] }
        flagstats = indcovflag_ch.map { tuple -> tuple[2] }
        combined_ch = metrics.mix(flagstats, FASTQC2_SE.out.logs).collect()
    }

    multiqc_config_ch = Channel.fromPath(params.multiqc_config, checkIfExists: true, type: 'file')

    MULTIQC2(combined_ch, multiqc_config_ch)
}
