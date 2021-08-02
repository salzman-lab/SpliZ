include { CONVERT_PARQUET } from '../../modules/local/convert_parquet' 

workflow PREPROCESS_TSV {
    take:
    ch_input

    main:
    CONVERT_PARQUET (
        ch_input,
        params.dataname
    )

    emit:
    pq = CONVERT_PARQUET.out.pq
}