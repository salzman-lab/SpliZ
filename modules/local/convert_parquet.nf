process CONVERT_PARQUET {
    tag "${dataname}"

    input:
    path tsv
    val dataname

    output:
    path "*.pq", emit: pq

    script:
    pq = "${tsv.baseName}.pq"
    """
    parquet_to_tsv.py \\
        --parquet ${pq} \\
        --tsv ${tsv} \\
        --reverse 
    """
}