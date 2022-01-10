process FIND_SPLIZ_SITES {
    tag "${params.dataname}"
    //label 'process_high_memory'
    publishDir "${params.outdir}/SpliZ_sites",  
        mode: "copy", 
        pattern: "*.tsv"
    
    label 'process_medium'

    input:
    path perm_pvals
    val libraryType
    path geneMat_samplesheet

    output:
    path first_evec     , emit: first_evec
    path second_evec    , emit: second_evec
    path third_evec     , emit: third_evec

    script:
    param_stem          = perm_pvals.baseName

    first_evec          = "first_evec_${param_stem}.tsv"
    second_evec         = "second_evec_${param_stem}.tsv"
    third_evec          = "third_evec_${param_stem}.tsv"

    """
    find_SpliZ_sites.R \\
        ${perm_pvals} \\
        ${first_evec} \\
        ${second_evec} \\
        ${third_evec} \\
        ${libraryType} \\
        ${geneMat_samplesheet}

    """

}
