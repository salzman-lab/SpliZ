#!/usr/bin/env python

import argparse
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description="convert parquet to tsv")
  parser.add_argument("--tsv",help="name to save tsv")
  parser.add_argument("--dataname",help="Dataname/basename of the input file")
  args = parser.parse_args()
  return args

def main():
  args = get_args()
  full_df = pd.read_csv(args.tsv, sep = "\t")

  df = full_df[full_df['called'] == True]

  for i, x in df.groupby('chrR1A'):
    outname = "{}_{}.pq".format(i, args.dataname)
    x.to_parquet(outname)

  
main()