process SUMMARIZE_RESULTS {
    tag "${params.dataname}"
    label 'process_low'
    publishDir "${params.outdir}",  
        mode: "copy", 
        pattern: "*.tsv"
    publishDir "${params.outdir}/logs", 
        mode: 'copy', 
        pattern: '*.log'

    input:
    path perm_pvals
    path first_evec
    path second_evec
    path third_evec
    val splizvd
    val group_col
    val sub_col

    output:
    path outname        , emit: summary
    path outname_log    , emit: log

    script:
    dataname            = splizvd[0]
    param_stem          = splizvd[1]
    splizvd_tsv         = splizvd[2]

    outname             = "summary_${dataname}_${group_col}-${sub_col}_${param_stem}.tsv"
    outname_log         = "summarize_results.log"

    """
    final_summary.py \\
        --perm_pvals ${perm_pvals} \\
        --first_evec ${first_evec} \\
        --second_evec ${second_evec} \\
        --third_evec ${third_evec} \\
        --splizvd ${splizvd_tsv} \\
        --group_col ${group_col} \\
        --sub_col ${sub_col} \\
        --outname ${outname} \\
        --outname_log ${outname_log}
    """
}