process ANN_SPLICES {
    tag "${params.dataname}"

    memory '800 GB'
    time '1h'

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