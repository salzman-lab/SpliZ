#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path
import pysam

def get_args():
  parser = argparse.ArgumentParser(description="merge class input files")
  parser.add_argument("--input_file", help="Metadata file")
  parser.add_argument("--meta", help="Metadata file")
  parser.add_argument("--outname", help="Output file name")
  parser.add_argument("--bam", help="bam file to extract intron reads from")

  parser.add_argument("--libraryType")

  args = parser.parse_args()
  return args

def get_intron_ret_df(bed_res, bam_file):
  df_bed = pd.read_csv(bed_res,sep="\t",names=["bed_chr","bed_start","bed_end","bed_name","bed_strand","bam_chr","bam_start","bam_end","read_id","mapq","bam_strand"])
  df_bed["len"] = df_bed["bam_end"] - df_bed["bam_start"]
  df_bed["num_readid"] = df_bed["read_id"].map(df_bed["read_id"].value_counts())
  read_len = df_bed.drop_duplicates("read_id")["len"].value_counts().index[0]
  df_bed = df_bed[(df_bed["num_readid"] == 1) & (df_bed["len"] == read_len)]
  read_dict = {k : v for k, v in zip(df_bed["read_id"],df_bed["bed_name"])}
  alignFile = pysam.AlignmentFile(bam_file)
  out = {'id' : [], 'refName_ABR1' : [], 'UMI' : [], 'barcode' : [], 'primaryR1' : [], 'read_strandR1' : [], 'juncPosR1A' : [], 'juncPosR1B' : [], 'geneR1A' : [], 'geneR1B' : [], 'chrR1A' : [], 'chrR1B' : []}
  
  # parse out "intron junctions"
  strand_dict={True : "-", False : "+"}
  for bam_read in (alignFile.fetch(until_eof=True)):
    if bam_read.query_name in read_dict.keys():
      seqname = alignFile.get_reference_name(bam_read.tid)
      refName = read_dict[bam_read.query_name]
      out["refName_ABR1"].append(refName)
      out["id"].append(bam_read.query_name)
      out["primaryR1"].append(not bam_read.is_secondary)
      out["juncPosR1A"].append(int(refName.split(":")[2]))
      out["juncPosR1B"].append(int(refName.split(":")[5]))
      out["chrR1A"].append(refName.split(":")[0])
      out["chrR1B"].append(refName.split(":")[3].split("|")[-1])
      out["geneR1A"].append(refName.split(":")[1])
      out["geneR1B"].append(refName.split(":")[4])
      out["read_strandR1"].append(strand_dict[bam_read.is_reverse])
      try:
        out["barcode"].append(bam_read.get_tag("CB"))
      except:
        out["barcode"].append(bam_read.get_tag("CR"))
      try:
        out["UMI"].append(bam_read.get_tag("UB"))
      except:
        out["UMI"].append(bam_read.get_tag("UR"))
  return pd.DataFrame.from_dict(out)

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
    print("df shape",df.shape)

    # intron retention inclusion
    int_df = get_intron_ret_df(row["file"][:-11] + "txt", args.bam)
    print("intron df shape",int_df.shape)
    df["intron"] = False
    int_df["intron"] = True
    df = pd.concat([df,int_df])
    print(df.head())
    print("intron in cols","intron" in df.columns)

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
  print("intron in cols","intron" in full_df.columns)

  full_df["called"] = 1
  full_df["refName_newR1"] = full_df["refName_ABR1"]
  full_df.rename(columns={"geneR1A" : "geneR1A_uniq", "geneR1B" : "geneR1B_uniq"}, inplace=True)
  
  final_df = full_df[["refName_newR1","geneR1A_uniq","geneR1B_uniq", "juncPosR1A","juncPosR1B","chrR1A","chrR1B","numReads","cell_id","intron"]]
  print("intron in cols","intron" in final_df.columns)
 
  meta = pd.read_csv(args.meta, sep="\t") 
  final_df.drop([x for x in final_df.columns if x in meta.columns and x != "cell_id"], inplace=True, axis=1)
  print("intron in cols","intron" in final_df.columns)

  merged = final_df.merge(meta, left_on="cell_id", right_on="cell_id", how = "left")
  print("intron in cols","intron" in merged.columns)

  merged.rename(columns={'cell_id': 'cell'}, inplace=True)
  merged.to_parquet(args.outname)


main()
