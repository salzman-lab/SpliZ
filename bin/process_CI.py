#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

def get_args():
  parser = argparse.ArgumentParser(description="merge class input files")
  parser.add_argument("--input_file", help="Metadata file")
  parser.add_argument("--meta", help="Metadata file")
  parser.add_argument("--outname", help="Output file name")
  parser.add_argument("--libraryType")

  args = parser.parse_args()
  return args

def main():
  args = get_args()

  file_list = pd.read_csv(args.input_file, header=None, names=['sample_ID','file'])

  file_list['sample_ID'] = file_list['sample_ID'].map(lambda x: x.lstrip('['))
  file_list['file'] = file_list['file'].map(lambda x: x.rstrip(']').lstrip(' '))

  dfs = []

  for index, row in file_list.iterrows():
    sample_ID = row['sample_ID']
    fn = Path(row['file'])
    df = pd.read_parquet(fn)

    # remove UMI duplicates by cell + junction
    df = df.drop_duplicates(["barcode","UMI","refName_ABR1"])
  
    df["barcode_refName"] = df["barcode"].astype(str) + df["refName_ABR1"]
  
    # count number of lines corresponding to the junction in the cell
    barcode_name_vc = df["barcode_refName"].value_counts()
    df["numReads"] = df["barcode_refName"].map(barcode_name_vc)

    # deduplicate by cell + junction
    df = df.drop_duplicates(["refName_ABR1","barcode"])

    # clean up barcode column
  
    if args.libraryType in ['10X',"SLS"]:
      df["barcode"] = df["barcode"].str.rstrip("-1")
      df["cell_id"] = sample_ID + "_" + df["barcode"].astype(str)
    elif args.libraryType == 'SS2':
      df['id'] = df['id'].str.split('.').str[0]
      df["cell_id"] = df["id"].astype(str)

  
    dfs.append(df)

  full_df = pd.concat(dfs)
  full_df["called"] = 1
  full_df["refName_newR1"] = full_df["refName_ABR1"]
  full_df.rename(columns={"geneR1A" : "geneR1A_uniq", "geneR1B" : "geneR1B_uniq"}, inplace=True)
  
  final_df = full_df[["refName_newR1","geneR1A_uniq","geneR1B_uniq", "juncPosR1A","juncPosR1B","chrR1A","chrR1B","numReads","cell_id"]]

  meta = pd.read_csv(args.meta, sep="\t") 
  final_df.drop([x for x in final_df.columns if x in meta.columns and x != "cell_id"], inplace=True, axis=1)

  merged = final_df.merge(meta, left_on="cell_id", right_on="cell_id", how = "left")

  merged.rename(columns={'cell_id': 'cell'}, inplace=True)
  merged.to_parquet(args.outname)


main()
