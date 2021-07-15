process CALC_SPLIZVD {
    tag "${dataname}"

    publishDir "${params.outdir}/SpliZ_values",  mode: "copy", pattern: "*.tsv"

    input:
    val rijk
    val svd_type

    output:
    tuple val(dataname), val(param_stem), path(outname_pq)  , emit: pq
    path "*.geneMat"                                        , emit: geneMats
    path outname_tsv                                        , emit: tsv

    script:
    dataname        = rijk[0]
    param_stem      = rijk[1]
    rijk_pq         = rijk[2]

    outname_pq      = "${dataname}_sym_SVD_${svd_type}_${param_stem}.pq"
    outname_tsv     = "${dataname}_sym_SVD_${svd_type}_${param_stem}_subcol.tsv"

    """
    svd_zscore.py \\
        --input ${rijk_pq} \\
        --svd_type ${svd_type} \\
        --outname_pq ${outname_pq} \\
        --outname_tsv ${outname_tsv} \\
    """
} 