#!/usr/bin/env python

import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import pickle
import pysam
import annotator
from light_utils import *
from tqdm import tqdm


def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--bams', nargs="+",required=True, help='bams to parse (either one or two for paired end)')
  parser.add_argument("--libraryType",help="Options: SS2, 10X")
  parser.add_argument("--annotator", required=True, help="the path to the annotator pickle file")
  parser.add_argument("--gtf", required=True, help="the path to the gtf file")
  parser.add_argument("--outname",help="Output file name")

  args = parser.parse_args()
  return args

def extract_info_align(cellranger, CI_dict, bam_read, suffix, bam_file, ann, UMI_bar, stranded_library, fill_char = np.nan, strand_dict={True : "-", False : "+"}):
  if UMI_bar:
    if cellranger:
      CI_dict["barcode"].append(bam_read.get_tag("CB"))
      CI_dict["UMI"].append(bam_read.get_tag("UR"))
    else:
      vals = bam_read.query_name.split("_")
      CI_dict["barcode"].append(vals[-2])
      CI_dict["UMI"].append(vals[-1])
  else:
    CI_dict["barcode"].append(fill_char)
    CI_dict["UMI"].append(fill_char)
  CI_dict["id"].append(bam_read.query_name)
  
  seqname = bam_file.get_reference_name(bam_read.tid)

  # if chromosome is numeric, prepend "chr"
  if str(seqname).isnumeric():
    seqname = "chr" + str(seqname)

  refName, chrA, geneA, posA, chrB, geneB, posB = readObj_refname(
    strand_dict[bam_read.is_reverse], 
    bam_read.cigarstring, 
    seqname, 
    bam_read.reference_start + 1, 
    ann, 
    fill_char, 
    stranded_library
  )
  CI_dict["refName_AB" + suffix].append(refName)
  CI_dict["chr{}A".format(suffix)].append(chrA)
  CI_dict["chr{}B".format(suffix)].append(chrB)
  CI_dict["gene{}A".format(suffix)].append(geneA)
  CI_dict["gene{}B".format(suffix)].append(geneB)
  CI_dict["juncPos{}A".format(suffix)].append(int(posA))
  if np.isnan(posB):
    CI_dict["juncPos{}B".format(suffix)].append(posB)
  else:
    CI_dict["juncPos{}B".format(suffix)].append(int(posB))
  strand_dict = {True : "-", False : "+"}
  CI_dict["read_strand{}".format(suffix)].append(strand_dict[bam_read.is_reverse])


  CI_dict["primary{}".format(suffix)].append(not bam_read.is_secondary)

  empty_cols = []
  for c in empty_cols:
    CI_dict[c].append(fill_char)
  return CI_dict
 
def extract_info_chim(CI_dict,bam_read1,bam_read2,suffix, bam_file, ann, UMI_bar, stranded_library, fill_char = np.nan):
  assert bam_read1.query_name == bam_read2.query_name
  sec_dict = {True: 0, False: 1}
  if UMI_bar:
    vals = bam_read1.query_name.split("_")
    CI_dict["barcode"].append(vals[-2])
    CI_dict["UMI"].append(vals[-1])
  else:
    CI_dict["barcode"].append(fill_char)
    CI_dict["UMI"].append(fill_char)
  reads = [bam_read1,bam_read2]
  halves = ["A","B"]
  CI_dict["id"].append(bam_read1.query_name)

  refName, chrA, geneA, posA, chrB, geneB, posB  = chim_refName([x.flag for x in reads], [x.cigarstring for x in reads], [x.reference_start + 1 for x in reads], [bam_file.get_reference_name(x.tid) for x in reads], ann, stranded_library)
  CI_dict["refName_AB" + suffix].append(refName)
  CI_dict["chr{}A".format(suffix)].append(chrA)
  CI_dict["chr{}B".format(suffix)].append(chrB)
  CI_dict["gene{}A".format(suffix)].append(geneA)
  CI_dict["gene{}B".format(suffix)].append(geneB)
  CI_dict["juncPos{}A".format(suffix)].append(int(posA))
  CI_dict["juncPos{}B".format(suffix)].append(int(posB))
  for i in range(2):

    CI_dict["primary{}{}".format(suffix,halves[i])].append(sec_dict[reads[i].is_secondary])
  return CI_dict


def get_final_df(cellranger, bam_files, j, suffixes, ann, UMI_bar, gtf, stranded_library):

  CI_dfs = []
  for i in range(len(bam_files)):
    if i == 1:
      read_ids = set(CI_dfs[0]["id"])
    else:
      read_ids = set()
    suffix = suffixes[i]
    col_bases = [ "juncPos", "gene", "chr"]
    columns = ["id", "refName_AB" + suffix, "UMI", "barcode", "primary" + suffix, "read_strand" + suffix]
    for c in col_bases:
      for l in ["A", "B"]:
        columns.append("{}{}{}".format(c,suffix,l))
    CI_dict = {c : [] for c in columns}
    count = 0
    first = False
    if i == 0:
      genomic_alignments = {}
    alignFile = pysam.AlignmentFile(bam_files[i])
    # columns
    #for bam_read in tqdm(alignFile.fetch(until_eof=True)):
    for bam_read in (alignFile.fetch(until_eof=True)):
      # require CB if this is cell ranger
      if ((not cellranger) | (bam_read.has_tag("CB"))):

        # make sure read is mapped
        if not bam_read.is_unmapped:
          if (i == 0) or (not bam_read.is_secondary and bam_read.query_name in read_ids):
            # it's a chimeric alignment and we need another line from it
            if bam_read.has_tag("ch") and not first:
                  prev_read = bam_read
                  first = True
            else:
  
              # add info from chimeric read
              if bam_read.has_tag("ch"):
                count += 1
  
                # note: removing chim for this test ONLY; uncomment after
                first = False
  
              # add info from align read
              elif "N" in bam_read.cigarstring:
                count += 1
                CI_dict = extract_info_align(cellranger, CI_dict, bam_read, suffix, alignFile, ann, UMI_bar, stranded_library)
  
              # save genomic alignment information
              else:
                if i == 0:
                  if bam_read.query_name not in genomic_alignments:
                    genomic_alignments[bam_read.query_name] = bam_read.get_tag("AS")
                  else:
                    genomic_alignments[bam_read.query_name] = max(bam_read.get_tag("AS"), genomic_alignments[bam_read.query_name])
                else:
                  CI_dict = extract_info_align(cellranger, CI_dict, bam_read, suffix, alignFile, ann, UMI_bar, stranded_library)

    CI_df = pd.DataFrame.from_dict(CI_dict)
    if i == 0:
      genomic_alignments = defaultdict(lambda: np.nan,genomic_alignments)

    CI_dfs.append(CI_df)
  if len(bam_files) == 2:
    final_df = pd.merge(left=CI_dfs[0],right=CI_dfs[1][[c for c in CI_dfs[1].columns if c not in ["UMI","barcode"]]],how="left",left_on="id",right_on="id")
    final_df["read_strand_compatible"] = 1
    final_df.loc[final_df["read_strandR1"] == final_df["read_strandR2"],"read_strand_compatible"] = 0
    final_df["location_compatible"] = final_df.apply(get_loc_flag,axis=1)
  else:
    final_df = CI_dfs[0]
  float_cols = ["primaryR1"]
  if len(bam_files) == 2:
    float_cols += ["juncPosR2A","juncPosR2B","primaryR2"]

  final_df = final_df[final_df["primaryR1"]]

  return final_df

def main():
  save = pysam.set_verbosity(0)
  
  args = get_args()

  bam_files = args.bams
  gtf = args.gtf

  annotator_path = args.annotator
  ann = pickle.load(open(annotator_path, "rb"))

  suffixes = ["R1","R2"]

  final_dfs = []
  
  n_rounds = len(bam_files)

  if args.libraryType == '10X':
    UMI_bar = True
    stranded_library = True
    cellranger = True

  elif args.libraryType == 'SS2':
    UMI_bar = False
    stranded_library = False
    cellranger = False

  for j in range(n_rounds):
    if j == 1:
      bam_files.reverse()
    primary = get_final_df(cellranger, bam_files, j, suffixes, ann, UMI_bar, gtf, stranded_library)
    final_dfs.append(primary)

  pd.concat(final_dfs, axis=0).reset_index(drop=True).to_parquet(args.outname)
  
  pysam.set_verbosity(save)



main()
