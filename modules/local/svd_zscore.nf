process SVD_ZSCORE {
    input:
    path rijk
    val temp
    val svd_type

    output:
    path "*.geneMat",   emit: geneMats
    path "*.svd",       emit: svd

    script:
    outname = "${rijk.baseName}.pq"
    """
    svd_zscore.py \\
        --rijk_file ${rijk} \\
        --outname ${outname} \\
        --temp ${temp} \\
        --svd_type ${svd_type} \\
    """
}