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
include { CONVERT_PARQUET       }   from   '../modules/local/convert_parquet' 
include { CALC_RIJK_ZSCORE      }   from   '../modules/local/calc_rijk_zscore'
include { CALC_SPLIZVD          }   from   '../modules/local/calc_splizvd'
include { PVAL_PERMUTATIONS     }   from   '../modules/local/pval_permutations'
include { FIND_SPLIZ_SITES      }   from   '../modules/local/find_spliz_sites'
include { SUMMARIZE_RESULTS     }   from   '../modules/local/summarize_results'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow SPLIZ {
    // Step 0a: Initialize input channel
    ch_pq = ch_input

    // Step 0b: If input file is not SICILIAN output, preprocess 

    // Step 0c: If input is tsv, convert to parquet
    if (convert_to_pq) {
        CONVERT_PARQUET (
            ch_input,
            params.dataname
        )
        ch_pq = CONVERT_PARQUET.out.pq
    } 

    // Step 1: Calculate RIJK zscore
    CALC_RIJK_ZSCORE (
        params.dataname,
        ch_pq,
        params.pin_S,
        params.pin_z,
        params.bounds,
        params.light,
        params.SICILIAN,
        params.grouping_level_2,
        params.grouping_level_1
    )
    
    // Step 2: Calculate SplizVD
    CALC_SPLIZVD (
        CALC_RIJK_ZSCORE.out.pq,
        params.svd_type      
    )
    
    if (!params.calc_SpliZ_only) {
        // Step 3: Calculate variance adjusted permutations
        PVAL_PERMUTATIONS (
            CALC_SPLIZVD.out.pq,
            params.n_perms,
            params.grouping_level_2,
            params.grouping_level_1
        )

        PVAL_PERMUTATIONS.out.perm_pvals
            .set{ ch_pval_permutations }

        // Step 4: Find SpliZ sites
        FIND_SPLIZ_SITES (
            CALC_SPLIZVD.out.geneMats, 
            ch_pval_permutations
        )

        // Step 5: Summarize results
        SUMMARIZE_RESULTS (
            ch_pval_permutations,
            FIND_SPLIZ_SITES.out.first_evec,
            FIND_SPLIZ_SITES.out.second_evec,
            FIND_SPLIZ_SITES.out.third_evec,
            CALC_SPLIZVD.out.tsv,
            params.grouping_level_2,
            params.grouping_level_1
        )
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/