include { CLASS_INPUT_10X       } from '../../modules/local/class_input_10X' 
include { CLASS_INPUT_SS2       } from '../../modules/local/class_input_SS2' 
include { PROCESS_CLASS_INPUT   } from '../../modules/local/process_class_input' 
include { ANN_SPLICES           } from '../../modules/local/ann_splices'

workflow CONVERT_BAM {
    take:
    ch_bam

    main:

    if ((params.libraryType == '10X') || (params.libraryType == "SLS")) {
        CLASS_INPUT_10X (
            ch_bam,
            params.dataname,
            params.libraryType,
            params.annotator_pickle,
            params.gtf
        )
        ch_light_class_input = CLASS_INPUT_10X.out.class_input
    } else if (params.libraryType == 'SS2') {
        CLASS_INPUT_SS2 (
                ch_bam,
                params.dataname,
                params.libraryType,
                params.annotator_pickle,
                params.gtf
            )
        ch_light_class_input = CLASS_INPUT_SS2.out.class_input
    }

    ch_light_class_input
        .collectFile(newLine: true) { files ->
            files.toString()
        }
        .set { ch_class_input }
    
    ch_class_input.view()
    
    PROCESS_CLASS_INPUT (
        ch_bam,
        ch_class_input,
        params.dataname,
        params.libraryType,
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
