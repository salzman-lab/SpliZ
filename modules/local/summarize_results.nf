process SUMMARIZE_RESULTS {
    tag "${params.dataname}"

    publishDir "${params.outdir}",  
        mode: "copy", 
        pattern: "*.tsv"
    publishDir "${params.outdir}/logs", 
        mode: 'copy', 
        pattern: '*.log'

    label 'process_low'

    input:
    path perm_pvals
    val param_stem
    val dataname
    path first_evec
    path second_evec
    path third_evec
    path splizvd_tsv
    val grouping_level_2
    val grouping_level_1

    output:
    path outname        , emit: summary
    path outname_log    , emit: log

    script:
    outname             = "summary_${dataname}_${grouping_level_2}-${grouping_level_1}_${param_stem}.tsv"
    outname_log         = "summarize_results.log"

    """
    final_summary.py \\
        --perm_pvals ${perm_pvals} \\
        --first_evec ${first_evec} \\
        --second_evec ${second_evec} \\
        --third_evec ${third_evec} \\
        --splizvd ${splizvd_tsv} \\
        --grouping_level_2 ${grouping_level_2} \\
        --grouping_level_1 ${grouping_level_1} \\
        --outname ${outname} \\
        --outname_log ${outname_log}
    """
}