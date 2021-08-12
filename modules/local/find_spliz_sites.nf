process FIND_SPLIZ_SITES {
    tag "${params.dataname}"
    //label 'process_high_memory'
    publishDir "${params.outdir}/SpliZ_sites",  
        mode: "copy", 
        pattern: "*.tsv"
    publishDir "${params.outdir}/logs", 
        mode: 'copy', 
        pattern: '*.log'
    
    memory '100 GB'
    time '1h'

    input:
    path ch_geneMats
    path perm_pvals

    output:
    path first_evec     , emit: first_evec
    path second_evec    , emit: second_evec
    path third_evec     , emit: third_evec
    path outname_log    , emit: log

    script:
    param_stem          = perm_pvals.baseName

    first_evec          = "first_evec_${param_stem}.tsv"
    second_evec         = "second_evec_${param_stem}.tsv"
    third_evec          = "third_evec_${param_stem}.tsv"
    outname_log         = "find_spliz_sites.log"

    """
    find_SpliZ_sites.R \\
        ${perm_pvals} \\
        ${first_evec} \\
        ${second_evec} \\
        ${third_evec} \\
        ${outname_log}
    """

}