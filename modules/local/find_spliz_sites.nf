process FIND_SPLIZ_SITES {
    tag "${dataname}"

    publishDir "${params.outdir}/SpliZ_sites",  mode: "copy", pattern: "*.txt"

    input:
    path ch_geneMats
    path svd_pvals

    output:
    path "tester.txt"

    script:
    """
    tester.R
    """

}