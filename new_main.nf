#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/spliz
========================================================================================
 nf-core/spliz Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/spliz
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CONVERT } from './modules/local/convert' 

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
        exit 1, "Incorrect input file type supplied, parquet file must be SICILIAN output."
    }
} else {
    exit 1, "Incorrect input file type supplied, accepted input formats are *.tsv and *.pq."
}


/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/
workflow {
    // If necessary, convert tsv to parquet.
    // Initialize input channel
    if (convert_to_pq) {
        CONVERT(input_file)
        ch_input = Channel.fromPath(CONVERT.out.pq)
    } else {
        ch_input = Channel.fromPath(input_file)
    }
}


