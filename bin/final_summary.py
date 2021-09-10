#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import logging

def get_args():
  parser = argparse.ArgumentParser(description="Create final summary file")
  parser.add_argument("--perm_pvals", help="Permutation pvalue file")
  parser.add_argument("--first_evec", help="First eigenvector file")
  parser.add_argument("--second_evec", help="Second eigenvector file")
  parser.add_argument("--third_evec", help="Third eigenvector file")
  parser.add_argument("--splizvd", help="SpliZVD file")
  parser.add_argument("--grouping_level_2", help="column to group the data by (e.g. ontology, compartment, tissue)", default="ontology")
  parser.add_argument("--grouping_level_1", help="subset data by this column before checking for differences (e.g. tissue, compartment)", default="dummy")
  parser.add_argument("--outname", help="Name of output file")
  parser.add_argument("--outname_log", help="Name of log file")
  parser.add_argument("--numGenes")
  parser.add_argument("--dataname")
  args = parser.parse_args()
  return args


def main():
  args = get_args()

  logging.basicConfig(
    filename = args.outname_log,
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

  logging.info("Starting")

  # load in data
  pval_df = pd.read_csv(args.perm_pvals, sep = "\t")
  
  splizsite_dfs = []
  evec_files = [args.first_evec, args.second_evec, args.third_evec]
  for evec_file in evec_files:
    splizsite_dfs.append(pd.read_csv(evec_file, sep="\t"))
  splizsite_df = pd.concat(splizsite_dfs,axis=0).drop_duplicates()
  
  df = pd.read_csv(args.splizvd, sep="\t")
  if (args.grouping_level_1 == "tiss_comp") & (args.grouping_level_1 not in df.columns):
    df["tiss_comp"] = df[args.grouping_level_1] + df[args.grouping_level_2]
  elif args.grouping_level_1 == "dummy":
    df["dummy"] = "dummy"

  # combine outputs
  out_dict = {"gene" : [],"grouping_level_1" : [], "grouping_level_2" : [],  "SpliZsites" : []}
  z_cols = ["scZ","svd_z0","svd_z1","svd_z2"]
  
  for z_col in z_cols:
    out_dict["{}_median".format(z_col)] = []
    out_dict["{}_pval".format(z_col)] = []
  
  for gene, gene_df in df.groupby("gene"):
    for tiss, tiss_df in gene_df.groupby(args.grouping_level_1):
      for ont, ont_df in tiss_df.groupby(args.grouping_level_2):
        out_dict["gene"].append(gene)
        out_dict["grouping_level_1"].append(tiss)
        out_dict["grouping_level_2"].append(ont)
        out_dict["SpliZsites"].append(",".join([str(x) for x in splizsite_df[splizsite_df["gene"] == gene]["end"]]))
        
        
        for z_col in z_cols:
  
          out_dict["{}_median".format(z_col)].append(ont_df[z_col].median())
          try:
            pval = pval_df[(pval_df["gene"] == gene) & ((pval_df["grouping_level_1"] == tiss) | (pval_df["grouping_level_1"].isna()))]["perm_pval_adj_{}".format(z_col)].iloc[0]
          except:
            pval = np.nan
          out_dict["{}_pval".format(z_col)].append(pval)
  out_df = pd.DataFrame.from_dict(out_dict)
  out_df = out_df.sort_values(["gene","grouping_level_1","scZ_median"])
  out_df.to_csv(args.outname, sep="\t", index=False)

  subset_df = out_df.dropna(subset=['SpliZsites']).sort_values(by='scZ_pval').reset_index()
  genes =  subset_df[['gene']].drop_duplicates().head(args.numGenes).reset_index()

  first_evec = pd.read_csv(args.first_evec, sep='\t')
  subset_first_evec = genes.merge(first_evec, on='gene')

  plotterFile_name = args.dataname + ".plotterFile"
  subset_first_evec.to_csv(plotterFile_name, index=False, sep='\t')

  logging.info("Completed")

main()