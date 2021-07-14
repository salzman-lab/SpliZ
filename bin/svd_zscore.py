#!/usr/bin python
  
import argparse
import datetime 
import numpy as np
import os
import pandas as pd
from scipy import linalg
import time
from tqdm import tqdm
import warnings

def get_args():
  parser = argparse.ArgumentParser(description="calculate splicing scores per gene/cell")
  parser.add_argument("--input", help="Name of the input file from rijk_zscore")
  parser.add_argument("--param_stem", help="Parameter string for the output file")
  parser.add_argument("--svd_type", choices=["normgene","normdonor"], help="Method of calculating matrix before SVD")
  args = parser.parse_args()
  return args

def main():
  args = get_args()

  df = pd.read_parquet(args.input)

  ##### PERFORM SVD ZSCORE CALCULATION #####

  if args.svd_type == "normgene":
    zcontrib_col = "zcontrib"
  elif args.svd_type == "normdonor":
    
    for let in ["A","B"]:
      # find number of reads per donor (or acceptor) per cell
      df["cell_gene_pos" + let] = df["cell_gene"] + df["juncPosR1" + let].astype(str)
      df["n.g_pos" + let] = df.groupby("cell_gene_pos" + let)["numReads"].transform("sum")
      
      # normalize on a donor/acceptor rather than a gene basis
      # TRY OUT NOT SQRT-ING denominator as normalization
      df["zcontrib_posnorm" + let] = df["numReads"] * df["nSijk" + let] / df["n.g_pos" + let]
    
    zcontrib_col = "zcontrib_posnorm"

  letters = ["A","B"]
  for let in letters:

    # replace NANs with zeros
    df["zcontrib{}_rep".format(let)] = df[zcontrib_col + let].fillna(0)

    # create label for each junction + donor/acceptor
    df["str_juncPosR1" + let] = df["juncPosR1" + let].astype(int).astype(str) + "_" + let
    df["cell_gene_pos" + let] = df["cell"] + df["geneR1A_uniq"] + df["juncPosR1" + let].astype(str)

    # get sum of zcontribs for the given cell and splice site
    df["summed_zcontrib" + let] = df.groupby("cell_gene_pos" + let)["zcontrib{}_rep".format(let)].transform('sum')

  k = 3 # number of components to include
  loads = {"f{}".format(i) : {} for i in range(k)}
  zs = {"svd_z{}".format(i) : {} for i in range(k)}
  
  for gene, gene_df in tqdm(df.groupby("geneR1A_uniq")):
    
    # get zcontrib matrix
    gene_mats = []
    for let in letters:
      gene_mat = gene_df.drop_duplicates("cell_gene_pos" + let).pivot_table(index="cell_gene",columns="str_juncPosR1{}".format(let),values="summed_zcontrib" + let,fill_value=0)

      gene_mats.append(gene_mat)
    gene_mat = gene_mats[0].merge(gene_mats[1],on="cell_gene")

    # mean-normalize the rows
    gene_mat = gene_mat.subtract(gene_mat.mean(axis=1),axis=0)
    
    # calculate svd
    u, s, vh = linalg.svd(gene_mat,check_finite=False,full_matrices=False)
    
    if len(s) >= k:
      # calculate new z scores based on svd
      new_zs = gene_mat.dot(np.transpose(vh[:k,:]))
  
      # calculate load on each component
      load = np.square(s)/sum(np.square(s))
  
      # save new zs and fs in dictionaries to save later
      for i in range(k):
        loads["f{}".format(i)][gene] = load[i]
        zs["svd_z{}".format(i)].update(pd.Series(new_zs[i].values,index=new_zs.index).to_dict())
  
      # save loadings
      v_out = pd.DataFrame(vh,columns=gene_mat.columns)
      gene_mat_name = "{}_{}_S_{}.geneMat".format(gene, args.dataname, args.param_stem)
      v_out.to_csv(gene_mat_name, index=False, sep = "\t")
      
  for i in range(k):
    df["f{}".format(i)] = df["geneR1A_uniq"].map(loads["f{}".format(i)])
    df["svd_z{}".format(i)] = df["cell_gene"].map(zs["svd_z{}".format(i)])
  
  df["svd_z_sumsq"] = (df[["svd_z{}".format(i) for i in range(k)]]**2).sum(axis=1)

  sub_cols = ["cell","geneR1A_uniq","scZ","svd_z_sumsq","n.g_A","n.g_B"] + ["f{}".format(i) for i in range(k)] + ["svd_z{}".format(i) for i in range(k)] #+ velocity_cols
  if "ontology" in df.columns:
    sub_cols = sub_cols + ["tissue","compartment","free_annotation","ontology"]
  
  svd_out_name = "{}_sym_SVD_{}_{}.pq".format(args.dataname, args.svd_type, args.param_stem)
  df.to_parquet(svd_out_name)


main()
