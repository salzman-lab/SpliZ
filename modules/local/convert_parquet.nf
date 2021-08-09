process CONVERT_PARQUET {
    tag "${params.dataname}"

    input:
    path tsv

    output:
    path "*.pq",    emit: pq

    script:
    pq = "${tsv.baseName}.pq"
    """
    parquet_to_tsv.py \\
        --parquet ${pq} \\
        --tsv ${tsv} \\
        --reverse 
    """
}