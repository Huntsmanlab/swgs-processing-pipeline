process WCX_REF_CREATE {
    tag "Creating WisecondorX reference"

    container 'docker://sofvdvel/wisecondorx:0.1'	
    cpus 16
    memory '32 G'
    time '12h'

	input:
    path(normals)
    val(binsize)
    path(bams)

    output:
    tuple val(binsize), path("wx_references")

    script: 
    """
        mkdir -p wx_references

        if [ ${params.wcx_ref_from_bams} = true ]; then
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
