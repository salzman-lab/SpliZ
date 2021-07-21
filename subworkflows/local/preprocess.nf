include { CONVERT_PARQUET } from '../../modules/local/convert_parquet' 

workflow PREPROCESS {
    take:
    ch_input
    convert_to_pq

    main:
    // Step 1: Initialize input channel
    ch_pq = ch_input

    // Step 2: If input file is not SICILIAN output, preprocess 

    // Step 3: If input is tsv, convert to parquet
    if (convert_to_pq) {
        CONVERT_PARQUET (
            ch_input,
            params.dataname
        )
        ch_pq = CONVERT_PARQUET.out.pq
    } 

    emit:
    ch_pq = ch_pq
}