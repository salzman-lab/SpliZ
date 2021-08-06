include { CLASS_INPUT           } from '../../modules/local/class_input' 
include { PROCESS_CLASS_INPUT   } from '../../modules/local/process_class_input' 
include { ANN_SPLICES           } from '../../modules/local/ann_splices'

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

    emit:
    tsv = ANN_SPLICES.out.tsv
}