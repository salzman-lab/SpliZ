include { CALC_SPLIZVD }   from   '../../modules/local/calc_splizvd'

workflow SPLIZ {
    take:
    ch_input

    main:

    def suff_light      = params.light     ? "_light" : ""
    def suff_SICILIAN   = params.SICILIAN  ? "_SICILIAN" : ""
    def suff_rank_quant = params.rank_quant == 0 ? "" : "_r_${params.rank_quant}"
    
    def isLight         = params.light     ? "1" : "0"
    def isSICILIAN      = params.SICILIAN  ? "1" : "0"

    param_stem          = "S_${params.pin_S}_z_${params.pin_z}_b_${params.bounds}_r_${suff_rank_quant}${suff_light}${suff_SICILIAN}"

    // Step 1: Calculate RIJK zscore
    CALC_SPLIZVD (
        ch_input,
        param_stem,
        params.dataname,
        params.pin_S,
        params.pin_z,
        params.bounds,
        params.svd_type,
        params.grouping_level_1,
        params.grouping_level_2,
        isLight,
        isSICILIAN,
        params.rank_quant
    )

    emit:
    geneMat_samplesheet = CALC_SPLIZVD.out.matSheet
    splizvd_tsv         = CALC_SPLIZVD.out.tsv
    splizvd_pq          = CALC_SPLIZVD.out.pq
    param_stem          = param_stem
}
