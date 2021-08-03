include { CALC_RIJK_ZSCORE      }   from   '../../modules/local/calc_rijk_zscore'
include { CALC_SPLIZVD          }   from   '../../modules/local/calc_splizvd'
include { MERGE_CHR             }   from   '../../modules/local/merge_chr'

workflow SPLIZ {
    take:
    ch_pq

    main:
    // Step 1: Calculate RIJK zscore
    CALC_RIJK_ZSCORE (
        ch_pq.flatten(),
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
        params.svd_type,
        params.grouping_level_2,
        params.grouping_level_1      
    )

    // Step 3: Re merge by chromosome
    CALC_SPLIZVD.out.tsv
        .collectFile(newLine: true) { files ->
            files.toString()
        }
        .set{ tsv_file_list }
    
    CALC_SPLIZVD.out.pq
        .collectFile(newLine: true) { files ->
            files.toString()
        }
        .set{ pq_file_list }

    MERGE_CHR (
        tsv_file_list,
        pq_file_list,
        params.dataname,
        CALC_SPLIZVD.out.param_stem,
        params.svd_type
    )

    emit:
    splizvd_geneMats    = CALC_SPLIZVD.out.geneMats
    splizvd_tsv         = MERGE_CHR.out.tsv
    splizvd_pq          = MERGE_CHR.out.pq
}