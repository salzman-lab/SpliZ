process ANN_SPLICES {
    tag "${params.dataname}"

    label 'process_medium'

    input:
    path pq
    path exon_pickle
    path splice_pickle

    output:
    path outname, emit: tsv

    script:
    outname = "${params.dataname}_ann_splices.tsv"
    """
    ann_splices.py \\
        --in_file ${pq} \\
        --out_file ${outname} \\
        --exon_pickle ${exon_pickle} \\
        --splice_pickle ${splice_pickle}
    """
}