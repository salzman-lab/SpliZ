include { CLASS_INPUT               } from '../../modules/local/class_input' 
include { PROCESS_CLASS_INPUT       } from '../../modules/local/process_class_input' 
include { ANN_SPLICES               } from '../../modules/local/ann_splices'
include { CONVERT_PARQUET           } from '../../modules/local/convert_parquet'
include { CONVERT_SPLIT_PARQUET     } from './../../modules/local/convert_split_parquet'

workflow CONVERT_BAM {
    take:
    ch_bam

    main:

    CLASS_INPUT (
        ch_bam,
        params.dataname,
        params.libraryType,
        params.annotator_pickle,
        params.gtf
    )

    CLASS_INPUT.out.class_input
        .collectFile(newLine: true) { files ->
            files.toString()
        }
        .view()
        .set { ch_class_input }
    
    PROCESS_CLASS_INPUT (
        ch_class_input,
        params.dataname,
        params.meta
    )

    ANN_SPLICES (
        PROCESS_CLASS_INPUT.out.pq,
        params.exon_pickle,
        params.splice_pickle
    )

    CONVERT_SPLIT_PARQUET (
        ANN_SPLICES.out.tsv
    )

    emit:
    pq = CONVERT_SPLIT_PARQUET.out.pq
}