process CN_ICHOR {
    publishDir params.output_directory, mode:'copy', pattern: "relative_cn/ichor/**"
    tag "Completing copy number analysis via ichor"
	errorStrategy 'retry', maxRetries: 5
    container 'docker://huntsmanlab/ichorcna:latest' 
    cpus 8
    memory '16 G'
    time '20h'
    
	input: 
    tuple val(sample_id), path(bam), path(bai), path(bam_txt), val(bam_type), val(bin_size)

    output:
    tuple path("relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.correctedDepth.txt"), path("relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.seg"), val(sample_id), val(bam_type), val(bin_size), emit: results
    path "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.cna.seg"
    path "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.correctedDepth.txt"
    path "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.seg"
    path "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.seg.txt"
    path "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/**/*"
 
    script:
    def bin_size_bases = bin_size.toInteger() * 1000
    def centromere_file
    if ("${params.reference_genome_version}" == "hg19" || "${params.reference_genome_version}" == "GRCh37") {
    centromere_file = "GRCh37.p13_centromere_UCSC-gapTable.txt"
    } else {
    centromere_file = "GRCh38.GCA_000001405.2_centromere_acen.txt"
    } 
    """
    mkdir -p "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/"
	bam_file=`readlink -f "${bam}"`	
	bamindex_file=`readlink -f "${bai}"`
	flock --verbose --wait 30 -E 3 "\${bamindex_file}.lock" cp -n \${bamindex_file} \${bamindex_file%.bai}.bam.bai    

    readCounter --window "${bin_size_bases}" --quality 20 \
    --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" \
    \${bam_file} > "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.readcounts.wig"
    
    Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R \
    --id "${sample_id}" \
    --WIG "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}/${sample_id}.readcounts.wig" \
    --gcWig "/usr/local/bin/ichorCNA/inst/extdata/gc_${params.reference_genome_version}_${bin_size}kb.wig" \
    --mapWig "/usr/local/bin/ichorCNA/inst/extdata/map_${params.reference_genome_version}_${bin_size}kb.wig" \
    --centromere "/usr/local/bin/ichorCNA/inst/extdata/${centromere_file}" \
    --ploidy "${params.ichor_ploidy}"  --normal "${params.ichor_normal}" \
    --maxCN "${params.ichor_maxCN}" \
    --minMapScore 0.75 \
    --includeHOMD True \
    --estimateScPrevalence False --txnE 0.9999 --txnStrength 10000 \
    --fracReadsInChrYForMale 0.002 \
    --outDir "relative_cn/ichor/${bam_type}/${bin_size}kb/${sample_id}"
    """
}
