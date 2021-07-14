process PERMUTATIONS {
    tag "${dataname}"

    publishDir "${params.outdir}",  mode: "copy", pattern: "*.tsv"

    input:
    

    output:


    script:
   

    """
    svd_zscore.py \\
        --input ${rijk_file} \\
        --dataname ${dataname} \\
        --param_stem ${param_stem} \\
        --svd_type ${svd_type} \\
    """
} 