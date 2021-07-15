process CALC_RIJK_ZSCORE {
    tag "${dataname}"

    publishDir "${params.outdir}/SpliZ_values"  , mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/logs"          , mode: 'copy', pattern: '*.log'

    input:
    val dataname
    path pq 
    val pin_S 
    val pin_z 
    val bounds 
    val light
    val SICILIAN

    output:
    tuple val(dataname), val(param_stem), path("*.pq")  , emit: pq
    path "*.tsv"                                        , emit: tsv    
    path "*.log"                                        , emit: log                                    

    script:
    def suff_light      = light     ? "_light" : ""
    def suff_SICILIAN   = SICILIAN  ? "_SICILIAN" : ""
    
    def isLight         = light     ? "0" : "1"
    def isSICILIAN      = SICILIAN  ? "0" : "1"

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
        --outname_pq ${outname_pq} \\
        --outname_tsv ${outname_tsv} \\
        --outname_log ${outname_log} \\
    """
}