process SVD_ZSCORE {
    publishDir("${params.resultsDir}/SVD_${svd_type}", mode: copy)

    input:
    path rijk
    val svd_type

    output:
    path "*.geneMat",   emit: geneMats
    path "*.pq",        emit: pq

    script:
    dataname = rijk[0]
    param_stem = rijk[1]
    rijk_file = rijk[2]

    """
    svd_zscore.py \\
        --input ${rijk_file} \\
        --param_stem ${param_stem} \\
        --svd_type ${svd_type} \\
    """
} 