process CONVERT_PARQUET {
    tag "${params.dataname}"
    label 'process_low'

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