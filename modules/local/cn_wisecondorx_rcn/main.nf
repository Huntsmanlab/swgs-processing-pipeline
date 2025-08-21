process CN_WISECONDORX_RCN {
    tag "WisecondorX RCN analysis with bin size ${binsize}kb"
	publishDir params.output_directory, mode:'copy'

    container 'docker://sofvdvel/wisecondorx:0.1'	
    cpus 16
    memory '32 G'
    time '12h'

    input:
    tuple val(binsize), path(wx_ref), val(bam_type), val(bams)
	path(output_dir)

    output:
    path("relative_cn/wisecondorx/${bam_type}/${binsize}kb")
    val(binsize)
    val(bam_type)

    script:
    """
    mkdir -p "relative_cn/wisecondorx/${bam_type}/${binsize}kb"
    mkdir -p converts
    printf '%s\n' "${bams.join('\n')}" > bamfileslist.txt

    while read LINE; do

        SAMPLE=\$(basename "\${LINE%%.*}")

        WisecondorX convert \$LINE converts/\${SAMPLE}.${bam_type}.npz --binsize ${binsize}000

    done < bamfileslist.txt

    for FILE in converts/*; do

        SAMPLE=\$(basename "\${FILE%%.*}")

        mkdir -p relative_cn/wisecondorx/${bam_type}/${binsize}kb/\${SAMPLE}

        if [[ \$FILE == *"pe."* ]]; then
            WisecondorX predict \$FILE ${wx_ref}/*${binsize}kb*.pe.npz relative_cn/wisecondorx/${bam_type}/${binsize}kb/\${SAMPLE}/\${SAMPLE} --gender F --plot --bed
        else
            WisecondorX predict \$FILE ${wx_ref}/*${binsize}kb*.se.npz relative_cn/wisecondorx/${bam_type}/${binsize}kb/\${SAMPLE}/\${SAMPLE} --gender F --plot --bed
        fi


    done

    rm bamfileslist.txt
    """
}
