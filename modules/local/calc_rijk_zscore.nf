process CALC_RIJK_ZSCORE {
    tag "${params.dataname}"
    //label 'process_high_memory'
    publishDir "${params.outdir}/SpliZ_values", 
        mode:       'copy', 
        pattern:    '*.tsv'
    publishDir "${params.outdir}/SpliZ_values", 
        mode:       'copy', 
        pattern:    '*.pq'
    publishDir "${params.outdir}/logs", 
        mode:       'copy', 
        pattern:    '*.log'

    input:
    val dataname
    path pq 
    val pin_S 
    val pin_z 
    val bounds 
    val light
    val SICILIAN
    val grouping_level_2
    val grouping_level_1

    output:
    tuple val(dataname), val(param_stem), path("*.pq")  , emit: pq
    path "*.tsv"                                        , emit: tsv    
    path "*.log"                                        , emit: log                                    

    script:
    def suff_light      = light     ? "_light" : ""
    def suff_SICILIAN   = SICILIAN  ? "_SICILIAN" : ""
    
    def isLight         = light     ? "1" : "0"
    def isSICILIAN      = SICILIAN  ? "1" : "0"

    param_stem          = "S_${pin_S}_z_${pin_z}_b_${bounds}${suff_light}${suff_SICILIAN}"

    outname_pq          = "${dataname}_sym_${param_stem}.pq"
    outname_tsv         = "${dataname}_sym_${param_stem}_subcol.tsv"
    outname_log         = "calc_rijk_zscore.log"
    
    """
    rijk_zscore.py \\
        --parquet ${pq} \\
        --pinning_S ${pin_S} \\
        --pinning_z ${pin_z} \\
        --lower_bound ${bounds} \\
        --isLight ${isLight} \\
        --isSICILIAN ${isSICILIAN} \\
        --grouping_level_2 ${grouping_level_2} \\
        --grouping_level_1 ${grouping_level_1} \\
        --outname_pq ${outname_pq} \\
        --outname_tsv ${outname_tsv} \\
        --outname_log ${outname_log} 
    """
}