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
    
    param_stem = CALC_RIJK_ZSCORE.out.param_stem

    // Step 2: Re merge by chromosome
    CALC_RIJK_ZSCORE.out.pq
        .collectFile(newLine: true) { file ->
            file.toString()
        }
        .set{ ch_rijk_pq }
    
    // Step 3: Calculate SplizVD
    CALC_SPLIZVD (
        ch_rijk_pq,
        params.dataname,
        param_stem,
        params.svd_type,
        params.grouping_level_2,
        params.grouping_level_1      
    )

    emit:
    splizvd_geneMats    = CALC_SPLIZVD.out.geneMats
    splizvd_tsv         = CALC_SPLIZVD.out.tsv
    splizvd_pq          = CALC_SPLIZVD.out.pq
    param_stem          = param_stem
}