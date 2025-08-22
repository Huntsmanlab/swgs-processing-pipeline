process ICHOR_COMBINE {
    publishDir params.output_directory, mode:'copy'
    tag "Combining ichor results"
    cpus 1
    memory '8 G'
    time '2h'
    
	input:
    tuple val(bamtype), val(binsize), path(tsv_files)

    output:
    path "relative_cn/ichor/combined_results/ichor_combined_${binsize}kb_${bamtype}_rCN.tsv"

    script:
    """
    mkdir -p "relative_cn/ichor/combined_results/"
    head -n 1 ${tsv_files[0]} > "relative_cn/ichor/combined_results/ichor_combined_${binsize}kb_${bamtype}_rCN.tsv"
    for FILE in ${tsv_files}; do
	tail -n +2 \$FILE >> "relative_cn/ichor/combined_results/combined_${binsize}kb_${bamtype}_rCN.tsv"
    done
    """
}
