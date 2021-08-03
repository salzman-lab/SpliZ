process CONVERT_SPLIT_PARQUET {
    tag "${params.dataname}"
    //label 'process_high_memory'

    input:
    path tsv

    output:
    path "*.pq",    emit: pq

    script:
    dataname = tsv.baseName
    """
    convert_tsv_to_parquet.py \\
        --tsv ${tsv} \\
        --splitChr \\
        --basename ${basename}
    """
}