process PROCESS_CLASS_INPUT {

    publishDir "${params.outdir}/class_input", 
        mode: 'copy', 
        pattern: '*.pq'

    label 'process_medium'

    input:
    path class_input
    val dataname
    val libraryType
    path meta

    output:
    path "*.pq",    emit: pq

    script:
    outname = "${dataname}.pq"
    """
    process_CI.py \\
        --input_file ${class_input} \\
        --meta ${meta} \\
        --libraryType ${libraryType} \\
        --outname ${outname} 
    """
}