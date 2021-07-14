process RIJK_ZSCORE {
    input:
    val dataname
    path pq 
    val pin_S 
    val pin_z 
    val light 
    val unfilt
    val v2
    val bound 

    output:
    path "*.svd", emit: svd

    script:
    def suff_light  = light ? "_light" : ""
    def suff_unfilt = unfilt ? "_unfilt" : ""

    def light       = light ? "0" : "1"
    def unfilt      = unfilt ? "0" : "1"
    def v2          = v2 ? "0" : "1"

    outname = "${dataname}_sym_S_${pin_S}_z_${pin_z}_b_${bound}${suff_light}${suff_unfilt}.pq"
    """
    rijk_zscore.py \\
        --dataname ${dataname} \\
        --parquet ${pq} \\
        --pinning_S ${pin_S} \\
        --pinning_z ${pin_z} \\
        --light ${light} \\
        --unfilt ${unfilt} \\
        --v2 ${v2} \\
        --lower_bound ${bound} \\
        --outname ${outname}
    """
}