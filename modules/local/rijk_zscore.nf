process RIJK_ZSCORE {
    publishDir "${params.outdir}"       , mode: copy, pattern: '*.tsv'
    publishDir "${params.outdir}/logs"  , mode: copy, pattern: '*.log'

    input:
    val dataname
    path pq 
    val pin_S 
    val pin_z 
    val bound 
    val light
    val SICILIAN

    output:
    tuple (val(dataname), val(param_stem), path("*.pq")), emit: pq

    script:
    def suff_light  = light ? "_light" : ""
    def suff_unfilt = unfilt ? "_unfilt" : ""

    def isSICILIAN  = SICILIAN ? "0" : "1"

    param_stem = "${pin_S}_z_${pin_z}_b_${bound}${suff_light}${suff_unfilt}"
    outname = "${dataname}_sym_S_${param_stem}"
    """
    rijk_zscore.py \\
        --dataname ${dataname} \\
        --parquet ${pq} \\
        --pinning_S ${pin_S} \\
        --pinning_z ${pin_z} \\
        --lower_bound ${bound} \\
        --light ${light} \\
        --isSICILIAN ${isSICILIAN} \\
        --outname ${outname}
    """
}