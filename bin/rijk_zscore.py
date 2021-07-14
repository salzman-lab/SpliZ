#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
import logging
warnings.filterwarnings("ignore")

def get_args():
  parser = argparse.ArgumentParser(description="calculate splicing scores per gene/cell")
  parser.add_argument("--dataname", help="name of dataset to use")
  parser.add_argument("--parquet", help="input parquet file")
  parser.add_argument("--pinning_S", type=float, help="pinning level for S_ijks")
  parser.add_argument("--pinning_z", type=float, help="pinning level for zs")
  parser.add_argument("--lower_bound", type=int, help="only include cell/gene pairs the have more than this many junctional reads for the gene")
  parser.add_argument("--isLight", help="if included, don't calculate extra columns (saves time)")
  parser.add_argument("--isSICILIAN", help="Is SICILIAN input file")
  parser.add_argument("--outname", type=str, help="Name of the output file")
  args = parser.parse_args()
  return args

def prepare_df(df, let, rank_by_donor, rev_let, let_dict):
  
  # create donor identifier
  df["pos{}_group".format(let)] = df["junc{}".format(let)].astype(str) + df["gene"]
  df["rank_" + let_dict[let]] = df.groupby("pos{}_group".format(let))["junc{}".format(rev_let[let])].rank(method="dense")

  # remove consitutive splicing
  df["max_rank"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["rank_" + let_dict[let]].max())
  df = df[df["max_rank"] > 1]
  
  if not rank_by_donor:
    df["rank_" + let_dict[let]] = df.groupby("gene")["juncEnd"].rank(method="dense")
  
  return df

def calc_Sijk(df,let, pinning_S, let_dict):
  # calculate the average rank calculation per gene
  # same as this calculation (this one's slower): df["rank_mean"] = df.groupby("pos{}_group".format(let)).apply(lambda x: (x["numReads"] * x["rank_acc"])/x["numReads"].sum()).reset_index(level=0,drop=True)

  # number of reads with this donor across all cells
  df["sum_reads_group"] = df.groupby("pos{}_group".format(let))["numReads"].transform("sum")

  df["read_x_" + let_dict[let]] = df["numReads"] * df["rank_" + let_dict[let]]

  # the sum of acceptors for all reads in all cells with this donor
  df["num"] = df.groupby("pos{}_group".format(let))["read_x_" + let_dict[let]].transform("sum")

  # average acceptor for a read with this donor (donor has one value for this)
  df["rank_mean"]= df["num"] / df["sum_reads_group"]

  # sum squared difference in rank for ever read
  df["sq_diff"] = df["numReads"] * (df["rank_" + let_dict[let]] - df["rank_mean"])**2

  # Get the sum of these squared differences for each donor
  df["don_num"] = df.groupby("pos{}_group".format(let))["sq_diff"].transform("sum")

  # sum of squared differences normalized by total number of reads
  # changed to make it the sample standard deviation (added minus 1)
  df["don_sigma"] = df["don_num"] / (df["sum_reads_group"])

  # this is the S_ijk value (difference normalized by sd) - should be normal 0/1
  df["S_ijk_{}".format(let)] = (df["rank_" + let_dict[let]] - df["rank_mean"])/np.sqrt(df["don_sigma"])

  # round outlying S values
  low_quant = df["S_ijk_{}".format(let)].quantile(q=pinning_S)
  high_quant = df["S_ijk_{}".format(let)].quantile(q=1 - pinning_S)
  df["S_ijk_{}_unpinned".format(let)] = df["S_ijk_{}".format(let)]

  df.loc[df["S_ijk_{}".format(let)] < low_quant,"S_ijk_{}".format(let)] = low_quant
  df.loc[df["S_ijk_{}".format(let)] > high_quant,"S_ijk_{}".format(let)] = high_quant

  # correct for those with no variance
  df.loc[df["don_sigma"] == 0, "S_ijk_{}".format(let)] = 0
  df["n_sijk"] = df["numReads"]
  df.loc[df["don_sigma"] == 0,"n_sijk"] = 0
  
  return df

def normalize_Sijks(df,let):

  # calculate mean of SijkA's per gene
  df["n_s"] = df["numReads"] * df["S_ijk_" + let]
  df["num"] = df.groupby("gene")["n_s"].transform("sum")
  df["n_gene"] = df.groupby("gene")["numReads"].transform("sum")
  df["sijk{}_mean".format(let)] = df["num"] / df["n_gene"]

  # calculate standard deviation of SijkA's per gene
  df["sd_num"] = df["numReads"] * (df["S_ijk_" + let] - df["sijk{}_mean".format(let)])**2
  df["num"] = df.groupby("gene")["sd_num"].transform("sum")
  df["sijk{}_var".format(let)] = df["num"] / df["n_gene"]

  return df

def contains_required_cols(df):

  # Function to check if the input file contains the required columns for processing
  required_cols = ["juncPosR1A", "geneR1A_uniq", "juncPosR1B", "numReads", "cell", "splice_ann", "tissue", "compartment", "free_annotation", "refName_newR1", "called", "chrR1A"]
  df_cols = list(df.columns)
  if set(required_cols) == set(df_cols):
    return True
  else:
    return False

def main():
  args = get_args()
  light = bool(args.isLight)
  SICILIAN = bool(args.isSICILIAN)

  outname_log = "{}.log".format(args.outname)
  outname_tsv = "{}.tsv".format(args.outname)
  outname_pq = "{}.pq".format(args.outname)

  logging.basicConfig(
    filename = outname_log,
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')  
  
  logging.info("Starting")

  let_dict = {"Start" : "acc", "End" : "don"}

  logging.info("Begin reading in parquet")

  df = pd.read_parquet(
    args.parquet,
    columns=["juncPosR1A", "geneR1A_uniq", "juncPosR1B", "numReads", "cell", "splice_ann", "tissue", "compartment", "free_annotation", "refName_newR1", "called", "chrR1A"]
  )
  
  logging.info("Finished reading in parquet")

  logging.info("Input column check")

  if contains_required_cols(df): 
    logging.info("Passed input column check")
  else:
    logging.info("Failed input column check! Exiting")
    sys.exit()

  logging.info("Rename SICILIAN columns")

  cols_dict = {
    "geneR1A_uniq": "gene",
    "juncPosR1A": "juncStart",
    "juncPosR1B": "juncEnd"
  }
  df.rename(columns=cols_dict, inplace=True)

  if "missing_domains" in df.columns and not light:
    domain_breakdown = True
  else:
    domain_breakdown = False

  df.reset_index(drop=True,inplace=True)
  rank_by_donor = True

  if SICILIAN:
    df = df[df["called"] == 1]
  else:
    # only include junctions with more than 1 read in the dataset
    df["numReads_tot"] = df.groupby("refName_newR1")["numReads"].transform("sum")
    df = df[df["numReads_tot"] > 1]

  # use second location gene name if first is unknown

  df["geneR1B_uniq"] = df["refName_newR1"].str.split("|").str[1].str.split(":").str[1]
  idx = df[(df["gene"].isin(["unknown",""])) | (df["gene"].isna())].index
  df.loc[idx,"gene"] = df.loc[idx,"geneR1B_uniq"]

  bin_size = 100000
  # bin unknown genes
  idx = df[(df["gene"] == "") | (df["gene"] == "unknown") | (df["gene"].isna())].index
  df.loc[idx,"gene"] = "unknown_" + df["chrR1A"].astype(str) + "_" + (df.loc[idx]["juncStart"] - df.loc[idx]["juncStart"] % bin_size).astype(str)
  print("replaced gene names",df[(df["gene"].isin(["unknown",""])) | (df["gene"].isna())].shape[0])

  logging.info("Replace with geneR1B")

  # get sign of gene to adjust z score
  sign_df = df.drop_duplicates("gene")
  sign_df["strandA"] = sign_df["refName_newR1"].str.split("|").str[0].str.split(":").str[3]
  sign_df["strandB"] = sign_df["refName_newR1"].str.split("|").str[1].str.split(":").str[3]
  idx = sign_df[sign_df["strandA"] == "?"].index
  sign_df.loc[idx,"strandA"] = sign_df.loc[idx,"strandB"]
  sign_df["sign"] = 1
  sign_df.loc[sign_df["strandA"] == "-","sign"] = -1
  sign_df[["gene","strandA","sign"]]
  sign_dict = pd.Series(sign_df.sign.values,index=sign_df.gene).to_dict()
  df["sign"] = df["gene"].map(sign_dict).fillna(1)

  logging.info("Get sign")

  df["cell_gene"] = df["cell"] + df["gene"]

  rev_let = {"Start" : "End", "End" : "Start"}

  if domain_breakdown:
    split_dict = {True : ["ann", "dom_ch"], False : ["unann", "dom_unch"]}
  else:
    split_dict = {True : ["ann"], False : ["unann"]}

  # remove constitutive splicing
  df["posA_group"] = df["juncStart"].astype(str) + df["gene"]
  df["posB_group"] = df["juncEnd"].astype(str) + df["gene"]

  df["rank_acc"] = df.groupby("posA_group")["juncEnd"].rank(method="dense")
  df["rank_don"] = df.groupby("posB_group")["juncStart"].rank(method="dense")

  df["max_rank_acc"] = df["posA_group"].map(df.groupby("posA_group")["rank_acc"].max())
  df["max_rank_don"] = df["posB_group"].map(df.groupby("posB_group")["rank_don"].max())

  # add domain columns
  letters = ["Start", "End"]
  for let in letters:

    if domain_breakdown:
      df["num_missing_" + let] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["missing_domains"].nunique())
      df["num_inserted_" + let] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["domain_insertions"].nunique())
      df["domain_changed_" + let] = (df["num_missing_" + let] + df["num_inserted_" + let]) > 0


  df = df[(df["max_rank_don"] > 1) | (df["max_rank_acc"] > 1)]

  logging.info("Remove constitutive")

  # require at least args.lower_bound nonconstitutive spliced reads
  df["noncon_count"] = df.groupby("cell_gene")["numReads"].transform("sum")
  df = df[df["noncon_count"] > args.lower_bound]

  full_df = df.copy()

  calc_dfs = {}

  for let in tqdm(letters):
    df = full_df
    # create donor identifier
    df = prepare_df(df, let, rank_by_donor, rev_let, let_dict)

    logging.info("Prepare df")
    df = calc_Sijk(df,let,args.pinning_S, let_dict)
    
    logging.info("Calculate Sijk")

    df = normalize_Sijks(df,let)

    logging.info("Normalize Sijk")
    
    # remove those with variance == 0
    df = df[df["sijk{}_var".format(let)] != 0]

    # calculate z score 
    df["n.g_" + let] = df.groupby("cell_gene")["numReads"].transform("sum")

    df["nSijk" + let] = (df["S_ijk_" + let] - df["sijk{}_mean".format(let)]) / np.sqrt(df["sijk{}_var".format(let)])
    df["mult"] = df["numReads"] * df["nSijk" + let]  / np.sqrt(df["n.g_" + let])
    df["z_" + let] = df["sign"] * df.groupby("cell_gene")["mult"].transform("sum")
    df["scaled_z_" + let] = df["z_" + let] / np.sqrt(df["n.g_" + let])

    logging.info("Calc z")

    ############## end modify Sijk ####################
    df["cell_gene_junc"] = df["cell_gene"] + df["refName_newR1"]

    if not light:
      # calculate the z score
      df["x_sijk"] = df["S_ijk_{}".format(let)] * df["n_sijk"]

      df["num"] = df.groupby("cell_gene")["x_sijk"].transform("sum")
      df["denom_sq"] = df.groupby("cell_gene")["n_sijk"].transform("sum")
  
      # get junction that "contributes the most" to the z score
      df["temp"] = df["x_sijk"] / np.sqrt(df["denom_sq"])
      df["temp_mag"] = abs(df["temp"])
      df["idxmax_z"] = df["cell_gene"].map(df.groupby("cell_gene")["temp_mag"].idxmax())
      map_df = df.loc[df["idxmax_z"],["cell_gene","refName_newR1","temp"]]
      df["junc_max_{}".format(let)] = df["cell_gene"].map(pd.Series(map_df.refName_newR1.values,index=map_df.cell_gene).to_dict())
      df["max_don_z_{}".format(let)] = df["cell_gene"].map(pd.Series(map_df.temp.values,index=map_df.cell_gene).to_dict())
  
    if args.pinning_z != 0:
      # round outlying z values
      low_quant = df["z_{}".format(let)].quantile(q=args.pinning_z)
      high_quant = df["z_{}".format(let)].quantile(q=1 - args.pinning_z)
      
      df.loc[df["z_{}".format(let)] < low_quant,"z_{}".format(let)] = low_quant
      df.loc[df["z_{}".format(let)] > high_quant,"z_{}".format(let)] = high_quant

    if not light:
      # break down z score by annotation
      for k,v in split_dict.items():
        df["num_{}".format(v[0])] = df["cell_gene"].map(df[df["splice_ann"] == k].groupby("cell_gene")["x_sijk"].sum())
  
        if domain_breakdown:
          df["num_{}".format(v[1])] = df["cell_gene"].map(df[df["domain_changed_" + let] == k].groupby("cell_gene")["x_sijk"].sum())
  
        for y in v:
  
          df["z_{}_{}".format(let,y)] = df["sign"] * df["num_{}".format(y)]/np.sqrt(df["denom_sq"])
      
          # round outlying z values
          low_quant = df["z_{}_{}".format(let,y)].quantile(q=args.pinning_z)
          high_quant = df["z_{}_{}".format(let,y)].quantile(q=1 - args.pinning_z)
      
          df.loc[df["z_{}_{}".format(let,y)] < low_quant,"z_{}_{}".format(let,y)] = low_quant
          df.loc[df["z_{}_{}".format(let,y)] > high_quant,"z_{}_{}".format(let,y)] = high_quant       

    calc_dfs[let] = df

  df = calc_dfs["Start"].merge(calc_dfs["End"],on="cell_gene_junc",how="outer",suffixes=("","_x"))
  
  logging.info("Merged")
  
  for cx in [x for x in df.columns if x.endswith("_x")]:
    c = cx[:-2]
    df.loc[df[c].isna(),c] = df.loc[df[c].isna(),cx]
  
  df.drop([x for x in df.columns if x.endswith("_x")],inplace=True,axis=1)

  # average two scores (negate one of them)

  grouped = df.groupby('gene')
  for let in letters:
    z_dict = pd.Series(calc_dfs[let]["z_" + let].values,index=calc_dfs[let].cell_gene).to_dict()
    df["z_" + let] = df["cell_gene"].map(z_dict)
    scz_dict = pd.Series(calc_dfs[let]["scaled_z_" + let].values,index=calc_dfs[let].cell_gene).to_dict()
    df["scaled_z_" + let] = df["cell_gene"].map(scz_dict)

  df["cov"] = df["gene"].map(grouped.apply(lambda x: x['z_Start'].cov(x['z_End'])))

  idx = df[df["z_Start"].isna()].index
  df.loc[idx,"z"] = -df.loc[idx,"z_End"]
  df.loc[idx,"scZ"] = -df.loc[idx,"scaled_z_End"]

  idx = df[df["z_End"].isna()].index
  df.loc[idx,"z"] = df.loc[idx,"z_Start"]
  df.loc[idx,"scZ"] = df.loc[idx,"scaled_z_Start"]

  idx = df[(~df["z_Start"].isna()) & (~df["z_End"].isna())].index
  df.loc[idx,"z"] = (df.loc[idx,"z_Start"] - df.loc[idx,"z_End"])/np.sqrt(2 )
  df.loc[idx,"scZ"] = (df.loc[idx,"scaled_z_Start"] - df.loc[idx,"scaled_z_End"])/np.sqrt(2 )

  logging.info("Avg z")
  
  if not light:
    # average two scores for split z
    for v in split_dict.values():
      for y in v:
        grouped = df.groupby('gene')
        df["cov_{}".format(y)] = df["gene"].map(grouped.apply(lambda x: x['z_Start_{}'.format(y)].cov(x['z_End_{}'.format(y)])))
  
        idx = df[df["z_Start_{}".format(y)].isna()].index
        df.loc[idx,"z_{}".format(y)] = -df.loc[idx,"z_End_{}".format(y)]
      
        idx = df[df["z_End_{}".format(y)].isna()].index
        df.loc[idx,"z_{}".format(y)] = df.loc[idx,"z_Start_{}".format(y)]
      
        idx = df[(~df["z_Start_{}".format(y)].isna()) & (~df["z_End_{}".format(y)].isna())].index
        df.loc[idx,"z_{}".format(y)] = (df.loc[idx,"z_Start_{}".format(y)] - df.loc[idx,"z_End_{}".format(y)])/np.sqrt(2) - df["cov_{}".format(y)]

  df["ontology"] = df["tissue"] + df["compartment"] + df["free_annotation"]
  
  df["n.g"] = df.groupby("cell_gene")["numReads"].transform("sum")
  df["scaled_z"] = df["z"] / np.sqrt(df["n.g"])

  for let in letters:
    df["zcontrib" + let] = df["numReads"] * df["nSijk" + let] / np.sqrt(df["n.g"])

  sub_cols = ["cell","gene","tissue","compartment","free_annotation","ontology","scZ","n.g_Start","n.g_End"] 

  df.drop_duplicates("cell_gene")[sub_cols].to_csv(outname_tsv, index=False,sep="\t")
  df.to_parquet(outname_pq)

  logging.info("Wrote files")

  logging.info("Completed")


main()