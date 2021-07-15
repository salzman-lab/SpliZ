process PVAL_PERMUTATIONS {
    tag "${params.dataname}"
    label 'process_medium'
    publishDir "${params.outdir}/variance_adjusted_permutations",  
        mode: "copy", 
        pattern: "*.tsv"

    input:
    val splizvd
    val n_perms
    val group_col
    val sub_col

    output:
    path outname_all_pvals,     emit: all_pvals
    path outname_perm_pvals,    emit: perm_pvals

    script:
    dataname            = splizvd[0]
    param_stem          = splizvd[1]
    splizvd_pq          = splizvd[2]

    outname_all_pvals   = "${dataname}_outdf_${group_col}-${sub_col}_${n_perms}_${param_stem}.tsv"
    outname_perm_pvals  = "${dataname}_pvals_${group_col}-${sub_col}_${n_perms}_${param_stem}.tsv"

    """
    variance_adjusted_permutations_bytiss.py \\
        --input ${splizvd_pq} \\
        --num_perms ${n_perms} \\
        --group_col ${group_col} \\
        --sub_col ${sub_col} \\
        --outname_all_pvals ${outname_all_pvals} \\
        --outname_perm_pvals ${outname_perm_pvals} \\
    """
} 