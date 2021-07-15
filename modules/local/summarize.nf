process SUMMARIZE {
    tag "${params.dataname}"

    publishDir "${params.outdir}",  mode: "copy", pattern: "*.tsv"

    input:
    path perm_pvals
    path first_evec
    path second_evec
    path third_evec
    path splizvd
    val group_col
    val sub_col

    output:
    path outname, emit: summary

    script:
    param_stem = perm_pvals.baseName

    outname = "summary_${param_stem}.tsv"

    """
    final_summary.py \\
        --perm_pvals ${perm_pvals} \\
        --first_evec ${first_evec} \\
        --second_evec ${second_evec} \\
        --third_evec ${third_evec} \\
        --splizvd ${splizvd} \\
        --group_col ${group_col} \\
        --sub_col ${sub_col} \\
        --outname ${outname} \\
    """
}