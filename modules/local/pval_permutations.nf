process PVAL_PERMUTATIONS {
    tag "${dataname}"

    publishDir "${params.outdir}/SpliZ_sites",  mode: "copy", pattern: "*.tsv"

    input:
    val svd
    val n_perms
    val group_col
    val sub_col

    output:
    path outname_all_pvals,     emit: all_pvals
    path outname_svd_pvals,     emit: svd_pvals

    script:
    dataname            = svd[0]
    param_stem          = svd[1]
    svd_pq              = svd[2]

    outname_all_pvals   = "${dataname}_outdf_${group_col}-${sub_col}_${n_perms}_${param_stem}.tsv"
    outname_svd_pvals   = "${dataname}_pvals_${group_col}-${sub_col}_${n_perms}_${param_stem}.tsv"

    """
    variance_adjusted_permutations_bytiss.py \\
        --input ${svd_pq} \\
        --num_perms ${n_perms} \\
        --group_col ${group_col} \\
        --sub_col ${sub_col} \\
        --outname_all_pvals ${outname_all_pvals} \\
        --outname_svd_pvals ${outname_svd_pvals} \\
    """
} 