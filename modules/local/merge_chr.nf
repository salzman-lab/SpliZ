process MERGE_CHR {
  tag "${params.dataname}"

publishDir "${params.outdir}/SpliZ_values",  
    mode: "copy"

  input:
  file tsv_file_list
  file pq_file_list
  val dataname
  val param_stem
  val svd_type

  output:
  path outname_tsv  , tsv
  path outname_pq   , pq

  script:
  outname_pq    = "${dataname}_sym_SVD_${svd_type}_${param_stem}.pq"
  outname_tsv   = "${dataname}_sym_SVD_${svd_type}_${param_stem}_subcol.tsv"
  """
  rm -f ${outname_tsv}
  rm -f ${outname_pq}

  cat ${tsv_file_list} |
    while read f; do
      cat \$f
    done >> ${outname_tsv}
  
  cat ${pq_file_list} |
    while read f; do
      cat \$f
    done >> ${outname_pq}
  """
}