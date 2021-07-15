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

// Initialise input channel
ch_input = Channel.fromPath(params.input_file)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CONVERT_PARQUET       }   from   '../modules/local/convert_parquet' 
include { CALC_RIJK_ZSCORE      }   from   '../modules/local/calc_rijk_zscore'
include { CALC_SPLIZVD          }   from   '../modules/local/calc_splizvd'
include { PVAL_PERMUTATIONS     }   from   '../modules/local/pval_permutations'
include { FIND_SPLIZ_SITES      }   from   '../modules/local/find_spliz_sites'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow SPLIZ {
    // Step 0: Initialize input channel
    ch_pq = ch_input
    if (convert_to_pq) {
        CONVERT_PARQUET (
            ch_input,
            params.dataname
        )
        ch_pq = CONVERT.out.pq
    } 

    // Step 1: Calculate RIJK zscore
    CALC_RIJK_ZSCORE (
        params.dataname,
        ch_pq,
        params.pin_S,
        params.pin_z,
        params.bounds,
        params.light,
        params.SICILIAN
    )
 
     // Step 2: Calculate SVD zscore
    CALC_SPLIZVD (
        RIJK_ZSCORE.out.pq,
        params.svd_type      
    )

    // Step 3: Calculate variance adjusted permutations
    PVAL_PERMUTATIONS (
        SVD_ZSCORE.out.pq,
        params.n_perms,
        params.group_col,
        params.sub_col
    )
}

/*
========================================================================================
    THE END
========================================================================================
*/