#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

def get_args():
  parser = argparse.ArgumentParser(description="calculate p values based on Romano method")
  parser.add_argument("--input", help="Name of the input file from svd_zscore")
  parser.add_argument("--num_perms", type=int,help="number of permutations to run for")
  parser.add_argument("--group_col", help="column to group the data by (e.g. ontology, compartment, tissue)", default="ontology")
  parser.add_argument("--sub_col", help="subset data by this column before checking for differences (e.g. tissue, compartment)", default="dummy")
  parser.add_argument("--outname_all_pvals", help="Name of output file")
  parser.add_argument("--outname_perm_pvals", help="Name of output File")
  parser.add_argument("--outname_log", help="Name of log File")
  args = parser.parse_args()
  return args


def calc_pval(var_df):
  
  # calculate the inner sum that's subtracted
  num = 0
  denom = 0
  for index, row in var_df.iterrows():
    num += row["num_cells_ont"]*row["ont_median"]/row["ont_var"]
    denom += row["num_cells_ont"]/row["ont_var"]
  const = num/denom

  # calculate the outer sum
  sum_vals = 0
  for index, row in var_df.iterrows():
    sum_vals += (row["num_cells_ont"]/row["ont_var"])*(row["ont_median"] - const)**2
    
  # return the chi^2 p value and the chi^2 statistic
  return 1 - stats.chi2.cdf(sum_vals , var_df.shape[0] - 1), sum_vals

def get_var_df(sub_df, z_col, adj_var, group_col):

  sub_df["num_cells_ont"] = sub_df[group_col].map(sub_df.groupby(group_col)["cell"].nunique())
  sub_df["ont_median"] = sub_df[group_col].map(sub_df.groupby(group_col)[z_col].median())
  sub_df["ont_var"] = sub_df[group_col].map(sub_df.groupby(group_col)[z_col].var())

  var_df = sub_df.drop_duplicates(group_col)[[group_col,"ont_median","num_cells_ont","ont_var"]]

  # don't need to remove cell types with variance 0 when we're adjusting variance
  if not adj_var:

    # remove ontologies with zero variance
    var_df = var_df[var_df["ont_var"] > 0]
  return var_df

def main():
  np.random.seed(123)
  alpha = 0.05

  args = get_args()

  logging.basicConfig(
    filename = args.outname_log,
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')  

  logging.info("Starting")

  #df_cols = ["gene", "cell", "scZ", "svd_z0", "svd_z1", "svd_z2", "cell_gene", "f0", "f1", "f2", "tissue", "compartment"]
  df_cols = ["gene", "cell", "scZ", "svd_z0", "svd_z1", "svd_z2", "cell_gene", "f0", "f1", "f2"]

  if args.sub_col.lower() != "dummy" and args.group_col not in df_cols:
    df_cols.append(args.group_col)  
    df_cols.append(args.group_col)

  if args.sub_col.lower == "dummy":
    df_cols.append(args.group_col)

  df = pd.read_parquet(
      args.input,
      columns=df_cols
  )
  df = df.drop_duplicates("cell_gene")
  df["tiss_comp"] = df["tissue"] + df["compartment"]

  # subset to ontologies with > 20 cells
  df["ontology_gene"] = df[args.group_col] + df["gene"]
  df["num_ont_gene"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell_gene"].nunique())
  df = df[df["num_ont_gene"] > 10]

  z_cols = ["scZ","svd_z0","svd_z1","svd_z2"]
  out = {"pval" : [], "gene" : [], "num_onts" : [],"z_col" : [],"max_abs_median" : [], "Tn1" : [], "sub_col" : []}
  
  var_adj = 0.1
  adj_var = True
  
  perm_pval = True
  
  if perm_pval:
    out["perm_pval"] = []
  
  df["dummy"] = "null"
  for tiss, tiss_df in df.groupby(args.sub_col):
    for gene, sub_df in tqdm(tiss_df.groupby("gene")):
      
      for z_col in z_cols:
    
        var_df = get_var_df(sub_df, z_col, adj_var, args.group_col)
        
        if var_df.shape[0] > 1:
          if adj_var:
            var_df["ont_var"] = var_df["ont_var"] + var_adj
          pval, Tn1 = calc_pval(var_df)
          out["pval"].append(pval)
          out["Tn1"].append(Tn1)
  
          out["gene"].append(gene)
          out["num_onts"].append(var_df.shape[0])
          out["z_col"].append(z_col)
          out["max_abs_median"].append((var_df["ont_median"].abs()).max())
          out["sub_col"].append(tiss)
          
          if perm_pval:
            sub_df_perm = sub_df.copy()
            if (pval < alpha):
              Tn1_dist = []
              # for i in range(args.num_perms):
              while len(Tn1_dist) < args.num_perms:
                sub_df_perm[args.group_col] = np.random.permutation(sub_df_perm[args.group_col])
                var_df = get_var_df(sub_df_perm, z_col, adj_var, args.group_col)
                if var_df.shape[0] > 1:
                  if adj_var:
                    var_df["ont_var"] = var_df["ont_var"] + var_adj
                  pval, Tn1_perm = calc_pval(var_df)
                  Tn1_dist.append(Tn1_perm)
              out["perm_pval"].append(len([x for x in Tn1_dist if x < Tn1])/args.num_perms)
            else:
              out["perm_pval"].append(np.nan)
  out_df = pd.DataFrame.from_dict(out)

  out_df["perm_pval_inv"] = 1 - out_df["perm_pval"] 
  out_df["perm_pval2"] = 2*out_df[["perm_pval","perm_pval_inv"]].min(axis=1)

  # adjust p values all together
  out_df["pval_adj"] = multipletests(out_df["pval"],alpha, method="fdr_bh")[1]
  out_df.loc[~out_df["perm_pval2"].isna(),"perm_pval2_adj"] = multipletests(out_df.loc[~out_df["perm_pval2"].isna(),"perm_pval2"], alpha, method = "fdr_bh")[1]

  out_df.to_csv(args.outname_all_pvals, sep="\t", index=False)

  out_df["gene_sub_col"] = out_df["gene"] + out_df["sub_col"]

  # reformat output
  new_out = {"gene" : [], "num_onts" : [], "sub_col" : []}
  for z_col in z_cols:
    new_out["chi2_pval_adj_" + z_col] = []
    new_out["perm_pval_adj_" + z_col] = []
    new_out["max_abs_median_" + z_col] = []
    new_out["perm_cdf_" + z_col] = []
  for gene_sub, gene_df in out_df.groupby("gene_sub_col"):
    new_out["gene"].append(gene_df["gene"].iloc[0])
    new_out["sub_col"].append(gene_df["sub_col"].iloc[0])
    new_out["num_onts"].append(gene_df["num_onts"].iloc[0])
    temp_z_cols = []
    for z_col, z_df in gene_df.groupby("z_col"):
      new_out["chi2_pval_adj_" + z_col].append(z_df["pval_adj"].iloc[0])
      new_out["perm_pval_adj_" + z_col].append(z_df["perm_pval2_adj"].iloc[0])
      new_out["max_abs_median_" + z_col].append(z_df["max_abs_median"].iloc[0])
      new_out["perm_cdf_" + z_col].append(z_df["perm_pval"].iloc[0])
      temp_z_cols.append(z_col)
    for z_col in [x for x in z_cols if x not in temp_z_cols]:
      new_out["chi2_pval_adj_" + z_col].append(np.nan)
      new_out["perm_pval_adj_" + z_col].append(np.nan)
      new_out["max_abs_median_" + z_col].append(np.nan)
      new_out["perm_cdf_" + z_col].append(np.nan) 
  new_out_df = pd.DataFrame.from_dict(new_out).sort_values("perm_pval_adj_scZ")

  # add frac from SVD for each gene
  df = df.drop_duplicates("gene")
  for i in range(3):
    frac_dict = pd.Series(df["f" + str(i)].values,index=df.gene).to_dict()
    new_out_df["f" + str(i)] = new_out_df["gene"].map(frac_dict)

  new_out_df.to_csv(args.outname_perm_pvals, sep="\t", index=False)
  
  logging.info("Completed")

main()