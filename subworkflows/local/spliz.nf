include { CALC_RIJK_ZSCORE      }   from   '../../modules/local/calc_rijk_zscore'
include { CALC_SPLIZVD          }   from   '../../modules/local/calc_splizvd'

workflow SPLIZ {
    take:
    ch_pq

    main:
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
        params.svd_type,
        params.grouping_level_2,
        params.grouping_level_1      
    )

    emit:
    splizvd_geneMats    = CALC_SPLIZVD.out.geneMats
    splizvd_tsv         = CALC_SPLIZVD.out.tsv
    splizvd_pq          = CALC_SPLIZVD.out.pq
}