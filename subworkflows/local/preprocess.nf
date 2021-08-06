include { CONVERT_PARQUET   } from './../../modules/local/convert_parquet'
include { CONVERT_BAM       } from './convert_bam'

workflow PREPROCESS {

    main:
    // Prepare 
    if (params.input_file) {

    }

    // Prepare inputs from SICILIAN output files
    if (params.SICILIAN) {   
        // Check that input_file is provided
        if (!params.input_file){
            exit 1, "No input_file provided."
        }

        // Stage input_file
        input_file = file(params.input_file)
        
        // Check if input file type is valid
        def is_valid_input_file = input_file.extension in ["tsv", "pq", "txt"]
        if (!is_valid_input_file) {
            exit 1, "Invalid input file type supplied, options are *.pq, *.txt, or *.tsv."
        } else { 
            // Initalize input channel
            ch_input = Channel.fromPath(params.input_file)
        }

    // Prepare inputs from non-SICILIAN bam files
    } else {
        // Initialize bam channel for bams specified in samplesheet
        if (params.samplesheet) {
            // Initialize bam channel for 10X bams specified in samplesheet
            ch_bam = Channel.fromPath(params.samplesheet)
                .splitCsv(header:false)
                .map { row ->
                    tuple( 
                        row[0],         // bam file sample_ID
                        file(row[1])    // bam file path 
                    )
                }   
        } else {
            exit 1, "No samplesheet provided."
        }

        // Check that bam channel is not empty
        ch_bam.ifEmpty (
            exit 1, "No bam files provied."
        )

        // Preprocess bam files for SpliZ pipeline
        CONVERT_BAM (
            ch_bam
        )

        // Initialize parquet channel for non-SICILIAN bam files
        ch_input = CONVERT_BAM.out.tsv
    }

    emit:
    input   = ch_input

}