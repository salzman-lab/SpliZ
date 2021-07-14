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
include { CONVERT       } from '../modules/local/convert' 
include { SVD           } from '../subworkflows/local/svd'
include { SPLIZ_SITES   } from '../subworkflows/local/spliz_sites'

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
    SVD (
        params.dataname,
        ch_pq,
        params.pin_S,
        params.pin_z, 
        params.bounds,
        params.light,
        params.SICILIAN,
        params.svd_type,
    )

    // Step 2: Identify SpliZ sites
    SPLIZ_SITES (
        SVD.out.geneMats,
    )

}

/*
========================================================================================
    THE END
========================================================================================
*/