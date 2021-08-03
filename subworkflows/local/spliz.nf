include { CALC_RIJK_ZSCORE      }   from   '../../modules/local/calc_rijk_zscore'
include { CALC_SPLIZVD          }   from   '../../modules/local/calc_splizvd'

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
        .flatten()
        .view()
        .map { file ->
            key = file.getName().toString().tokenize("_")[0]
            return tuple(key, file)
        }
        .groupTuple()
        .view()
        .collectFile { id, files ->
            [
                id,
                files.collect{ it.toString() }.join('\n') + '\n'
            ]
        }
        .set { channel_merge_list }
    
    CALC_SPLIZVD.out.geneMats
        .flatten()
        .map { file ->
            key = file.toString().tokenize("_")[0]
            return tuple(key, file)
        }
        .groupTuple()
        .view()
        .collectFile { id, files ->
            [
                id,
                files.collect{ it.toString() }.join('\n') + '\n'
            ]
        }
        .set { channel_merge_list }

    emit:
    splizvd_geneMats    = CALC_SPLIZVD.out.geneMats
    splizvd_tsv         = CALC_SPLIZVD.out.tsv
    splizvd_pq          = CALC_SPLIZVD.out.pq
}