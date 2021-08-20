include { CONVERT_PARQUET   } from './../../modules/local/convert_parquet'
include { CONVERT_BAM       } from './convert_bam'

workflow PREPROCESS {

    main:

    convert_bam = false

    if (params.input_file && params.samplesheet) {
        exit 1, "Invalid input, provide either input_file or samplesheet but not both."
    } else if (params.samplesheet) {
        if (params.SICILIAN) {
            exit 1, "Invalid input, SICILIAN inputs must be provided as input_file."
        } else {
            if (params.libraryType == '10X') {
                ch_bam = Channel.fromPath(params.samplesheet)
                    .splitCsv(header:false)
                    .map { row ->
                        tuple( 
                            row[0],         // bam file sample_ID
                            file(row[1])    // bam file path 
                        )
                    }
                convert_bam = true
            } else if (params.libraryType == 'SS2') {
                ch_bam = Channel.fromPath(params.samplesheet)
                    .splitCsv(header:false)
                    .map { row ->
                        tuple( 
                            row[0],         // bam file sample_ID
                            file(row[1]),   // R1 bam file path 
                            file(row[2])    // R2 bam file path    
                        )
                    }
                convert_bam = true
            }
        }
    } else if (params.input_file) {
        input_file = file(params.input_file)
        def is_valid_input_file = input_file.extension in ["tsv", "pq", "txt", "bam"]
        if (!is_valid_input_file) {
            exit 1, "Invalid input file type supplied, options are *.bam, *.pq, *.txt, or *.tsv."
        } 
        if (params.SICILIAN) {
            if (input_file.extension == "bam") {
                exit 1, "Invalid input, SICILIAN input must be a tsv, pq, or txt file."
            } else {
                ch_input = Channel.fromPath(params.input_file)
            }
        } else {
            if (input_file.extension == "bam") {
                if (!params.dataname) {
                    exit 1, "Must provide dataname for bam file."
                }
                ch_bam = Channel.fromPath(params.input_file)
                    .map { it ->
                        tuple(
                            params.dataname,
                            file(it)
                        )
                    }
                convert_bam = true
            } else {
                ch_input = Channel.fromPath(params.input_file)
            }
        }
    } else {
        exit 1, "No input_file or samplesheet provided."
    }
    
    if (convert_bam) {
        CONVERT_BAM (
            ch_bam
        )
        ch_input = CONVERT_BAM.out.tsv
    }

    emit:
    input = ch_input

}