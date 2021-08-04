process MERGE_CHR {
  tag "${params.dataname}"

publishDir "${params.outdir}/SpliZ_values",  
    mode: "copy"

  input:
  file tsv_file_list
  file pq_file_list
  val dataname
  val param_stem

  output:
  path "*.pq"   , emit: pq

  script:
  outname_pq    = "${dataname}_sym_${param_stem}.pq"
  outname_tsv   = "${dataname}_sym_${param_stem}_subcol.tsv"
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