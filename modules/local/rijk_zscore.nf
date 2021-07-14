process RIJK_ZSCORE {
    tag "${dataname}"

    publishDir "${params.outdir}"       , mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/logs"  , mode: 'copy', pattern: '*.log'

    input:
    val dataname
    path pq 
    val pin_S 
    val pin_z 
    val bound 
    val light
    val SICILIAN

    output:
    tuple val(dataname), val(param_stem), path("*.pq")  , emit: pq
    path "*.tsv"                                        , emit: tsv    
    path "*.log"                                        , emit: log                                    

    script:
    def suff_light      = light ? "light" : ""
    def suff_SICILIAN   = SICILIAN ? "SICILIAN" : ""
    
    def isLight         = light ? "0" : "1"
    def isSICILIAN      = SICILIAN ? "0" : "1"

    param_stem = "S_${pin_S}_z_${pin_z}_b_${bound}_${suff_light}_${suff_SICILIAN}"
    outname = "${dataname}_sym_${param_stem}"
    
    """
    rijk_zscore.py \\
        --dataname ${dataname} \\
        --parquet ${pq} \\
        --pinning_S ${pin_S} \\
        --pinning_z ${pin_z} \\
        --lower_bound ${bound} \\
        --isLight ${isLight} \\
        --isSICILIAN ${isSICILIAN} \\
        --outname ${outname}
    """
}