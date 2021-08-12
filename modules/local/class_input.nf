process CLASS_INPUT {
    tag "${params.dataname}"
    memory '400 GB'

    input:
    val bam
    val dataname
    val libraryType
    path annotator_pickle
    path gtf

    output:
    tuple val(sample_ID), path (outname),    emit: class_input

    script:
    if (libraryType == "10X") {
        sample_ID   = bam[0]
        bam_R1      = bam[1]
    } 
    if (libraryType == "SS2") {
        sample_ID   = bam[0]
        bam_R1      = bam[1]
        bam_R2      = bam[2]
    }

    outname         = "${sample_ID}.class_input"

    if (libraryType == "10X")
        """
        light_class_input_subcols.py \\
            --bams ${bam_R1} \\
            --libraryType ${libraryType} \\
            --annotator ${annotator_pickle} \\
            --gtf ${gtf} \\
            --outname ${outname} 
        """
    else if (libraryType == "SS2")
        """
        light_class_input_subcols.py \\
            --bams ${bam_R1} ${bam_R2} \\
            --libraryType ${libraryType} \\
            --annotator ${annotator_pickle} \\
            --gtf ${gtf} \\
            --outname ${outname} 
        """
}