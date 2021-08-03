process PROCESS_CLASS_INPUT {
    tag "${params.dataname}"
    label 'process_high'
    publishDir "${params.outdir}/class_input", 
        mode: 'copy', 
        pattern: '*.pq'

    input:
    path class_input
    val dataname
    path meta

    output:
    path "*.pq",    emit: pq

    script:
    outname = "${dataname}.pq"
    """
    process_CI.py \\
        --input_file ${class_input} \\
        --meta ${meta} \\
        --outname ${outname} 
    """
}