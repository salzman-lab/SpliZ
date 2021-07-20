/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
// Stage input file
input_file = file(params.input_file)

// Determine if file type needs to be converted to parquet.
if (input_file.extension == "tsv") {
    convert_to_pq = true
} else if (input_file.extension == "pq") {
    convert_to_pq = false
    if (!params.SICILIAN) {
        exit 1, "Invalid input file type supplied, parquet file must be SICILIAN output."
    }
} else {
    exit 1, "Invalid input file type supplied, accepted input formats are *.tsv and *.pq."
}

// Validate svd_type param
def is_valid_svd_type = params.svd_type in ["normgene", "normdonor"]
if (!is_valid_svd_type) {
    exit 1, "Invalid svd_type; options are 'normgene' and 'normdonor'."
}

// Initialise input channel
ch_input = Channel.fromPath(params.input_file)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { PREPROCESS        } from './../subworkflows/local/preprocess'
include { SPLIZ             } from './../subworkflows/local/spliz'
include { ANALYSIS          } from './../subworkflows/local/analysis'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow SPLIZ_PIPELINE {
    
    PREPROCESS (
        ch_input,
        convert_to_pq
    )

    SPLIZ (
        PREPROCESS.out.ch_pq
    )
    
    if (!params.calc_SpliZ_only) {
        ANLYSIS (
            SPLIZ.out.splizvd_geneMats,
            SPLIZ.out.splizvd_tsv,
            SPLIZ.out.splizvd_pq
        )
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/