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

ch_input = Channel.fromPath(params.input_file)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CONVERT } from '../modules/local/convert' 

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow SPLIZ {
    // If necessary, convert tsv to parquet.
    // Initialize input channel
    if (convert_to_pq) {
        CONVERT(ch_input)
        ch_pq = CONVERT.out.pq
    } else {
        ch_pq = ch_input
    }

    // Step 1: SVD calculation
    SVD(
        dataname, 
        ch_pq,
        pin_S,
        pin_Z,
        light,
        unfilt,
        v2,
        bound, 
        temp,
        svd_type
    )

}

/*
========================================================================================
    THE END
========================================================================================
*/