#!/usr/bin/env python

import argparse
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description="convert parquet to tsv")
  parser.add_argument("-p","--parquet",help="input parquet file")
  parser.add_argument("-o","--tsv",help="name to save tsv")
  parser.add_argument("--splitChr",action="store_true",help="split file by chromosome")
  args = parser.parse_args()
  return args

def main():
  args = get_args()

  df = pd.read_csv(args.tsv, sep = "\t")
  
  if args.splitChr:
    for i, x in df.groupby('chrR1A'):
      outfile = ("{}")
  else:
    df.to_parquet(args.parquet)

     
  
main()