process CLASS_INPUT_10X {
    tag "${params.dataname}"
    
    label 'process_high'

    input:
    tuple val(sample_ID), file(bam)
    val dataname
    val libraryType
    path annotator_pickle
    path gtf

    output:
    tuple val(sample_ID), path(outname),    emit: class_input

    script:
    outname = "${sample_ID}.class_input"

    """
    light_class_input_subcols.py \\
        --bams ${bam} \\
        --libraryType ${libraryType} \\
        --annotator ${annotator_pickle} \\
        --gtf ${gtf} \\
        --outname ${outname} 

    bedtools intersect -a ${sample_ID}.bed -b ${bam}  -wa -wb -f 1 >  ${sample_ID}.temp 
    awk '{if (\$10==255) print }' ${sample_ID}.temp > ${sample_ID}.txt
    """

}
