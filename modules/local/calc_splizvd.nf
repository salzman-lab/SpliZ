process CALC_SPLIZVD {
    tag "${params.dataname}"
    label 'process_high'
    publishDir "${params.outdir}/SpliZ_values",  
        mode: "copy", 
        pattern: "*.tsv"
    publishDir "${params.outdir}/SpliZ_values",  
        mode: "copy", 
        pattern: "*.pq"
    publishDir "${params.outdir}/logs", 
        mode: 'copy', 
        pattern: '*.log'
    
    input:
    val rijk
    val svd_type

    output:
    tuple val(dataname), val(param_stem), path(outname_pq)  , emit: pq
    tuple val(dataname), val(param_stem), path(outname_tsv) , emit: tsv                                 
    path "*.geneMat"                                        , emit: geneMats
    path "*.log"                                            , emit: log                                    

    script:
    dataname        = rijk[0]
    param_stem      = rijk[1]
    rijk_pq         = rijk[2]

    outname_pq      = "${dataname}_sym_SVD_${svd_type}_${param_stem}.pq"
    outname_tsv     = "${dataname}_sym_SVD_${svd_type}_${param_stem}_subcol.tsv"
    outname_log     = "calc_splizvd.log"

    """
    svd_zscore.py \\
        --input ${rijk_pq} \\
        --svd_type ${svd_type} \\
        --outname_pq ${outname_pq} \\
        --outname_tsv ${outname_tsv} \\
        --outname_log ${outname_log}
    """
} 