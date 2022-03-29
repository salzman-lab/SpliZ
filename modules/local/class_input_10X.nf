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

    samtools view -H ${bam} | grep -P "@SQ\tSN:"| sed 's/@SQ\tSN://'| sed 's/\tLN:/\t/' > genome.txt
    bedtools sort -g genome.txt -i ${sample_ID}.bed  > ${sample_ID}.sorted.bed
    bedtools intersect -sorted -g genome.txt -a ${sample_ID}.sorted.bed -b ${bam}  -wa -wb -f 1 >  ${sample_ID}.temp 
    awk '{if (\$10==255) print }' ${sample_ID}.temp > ${sample_ID}.txt
    """

}
