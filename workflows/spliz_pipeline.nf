/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check params with defined inputs
def is_valid_svd_type = params.svd_type in ["normgene", "normdonor"]
if (!is_valid_svd_type) {
    exit 1, "Invalid svd_type; options are 'normgene' and 'normdonor'."
}

def is_valid_libraryType = params.libraryType in ["SS2", "10X"]
if (!is_valid_libraryType) {
    exit 1, "Invalid libraryType; options are 'SS2' and '10X'."
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { PREPROCESS    } from './../subworkflows/local/preprocess'
include { SPLIZ         } from './../subworkflows/local/spliz'
include { ANALYSIS      } from './../subworkflows/local/analysis'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SPLIZ_PIPELINE {

    PREPROCESS ()

    
    SPLIZ (
        PREPROCESS.out.pq
    )
    
    if (params.run_analysis) {
        ANALYSIS (
            SPLIZ.out.param_stem,
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