process PVAL_PERMUTATIONS {
    tag "${params.dataname}"
    label 'process_medium'
    publishDir "${params.outdir}/variance_adjusted_permutations",  
        mode: "copy", 
        pattern: "*.tsv"
    publishDir "${params.outdir}/logs", 
        mode: 'copy', 
        pattern: '*.log'

    input:
    val splizvd
    val n_perms
    val grouping_level_2
    val grouping_level_1

    output:
    path outname_all_pvals      , emit: all_pvals
    path outname_perm_pvals     , emit: perm_pvals
    path outname_log            , emit: log

    script:
    dataname                    = splizvd[0]
    param_stem                  = splizvd[1]
    splizvd_pq                  = splizvd[2]

    outname_all_pvals           = "${dataname}_outdf_${grouping_level_2}-${grouping_level_1}_${n_perms}_${param_stem}.tsv"
    outname_perm_pvals          = "${dataname}_pvals_${grouping_level_2}-${grouping_level_1}_${n_perms}_${param_stem}.tsv"
    outname_log                 = "pval_permutations.log"

    """
    variance_adjusted_permutations_bytiss.py \\
        --input ${splizvd_pq} \\
        --num_perms ${n_perms} \\
        --grouping_level_2 ${grouping_level_2} \\
        --grouping_level_1 ${grouping_level_1} \\
        --outname_all_pvals ${outname_all_pvals} \\
        --outname_perm_pvals ${outname_perm_pvals} \\
        --outname_log ${outname_log}
    """
} 