process PLOT {
    tag "${params.dataname}"
    publishDir "${params.outdir}/plots/dotplots",  
        mode: "copy", 
        pattern: "*.dotplot"
    publishDir "${params.outdir}/plots/boxplots",  
        mode: "copy", 
        pattern: "*.boxplot"

    label 'process_medium'

    input:
    path plotterFile
    path splizvd_pq
    path domain
    path gtf
    val dataname
    val grouping_level_1
    val grouping_level_2

    output:
    path '*dotplot'     , emit: dotplot   
    path '*boxplot'     , emit: boxplot                                  
                               
    script:
    """
    plot.py \\
        --plotterFile ${plotterFile} \\
        --svd ${svd} \\
        --domain ${domain} \\
        --gtf ${gtf} \\
        --dataname ${dataname} \\
        --grouping_level_1 ${grouping_level_1} \\
        --grouping_level_2 ${grouping_level_2}
    """

} 