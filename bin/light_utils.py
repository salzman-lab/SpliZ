from collections import defaultdict
import math
import numpy as np
import pandas as pd
import re


def get_gene_id(row):
  if "gene_name" in row["attribute"]:
    return row["attribute"].split("gene_name")[-1].split('"')[1]
  elif ";gene=" in row["attribute"]:
    return row["attribute"].split(";gene=")[-1].split(";")[0]

def modify_refnames(CI, gtf_file, stranded_library):
  
  gtf_df = pd.read_csv(gtf_file,sep="\t",names=["seqname","source","feature","start","end","score","strand","frame","attribute"],comment="#")
  gtf_df["gene_name"] = gtf_df.apply(get_gene_id, axis=1)
  gtf_df = gtf_df[['seqname', 'strand','gene_name']]
  gene_strand_info = gtf_df.drop_duplicates().reset_index(drop=True)
  swap_names = False
  CI["HIR1B"] = CI["HIR1A"]
  CI_new = CI.drop_duplicates("refName_ABR1")
  

  CI_new["geneR1A"] = CI_new["geneR1A"].fillna("")
  CI_new["geneR1B"] = CI_new["geneR1B"].fillna("")
  CI_new.loc[CI_new["fileTypeR1"] == "Aligned","read_strandR1B"] = CI_new[CI_new["fileTypeR1"] == "Aligned"]["read_strandR1A"]
  CI_new["gene_strandR1A"] = CI_new["refName_ABR1"].str.split("|").str[0].str.split(":").str[-1]
  CI_new["gene_strandR1B"] = CI_new["refName_ABR1"].str.split("|").str[1].str.split(":").str[-1]
  CI_new["numgeneR1A"] = CI_new["geneR1A"].str.split(",").str.len()#.astype("Int32") # the number of overlapping genes on the R1A side
  CI_new[["numgeneR1A"]] = CI_new[["numgeneR1A"]].fillna(0)
  CI_new["numgeneR1B"] = CI_new["geneR1B"].str.split(",").str.len()#.astype("Int32") # the number of overlapping genes on the R1B side
  CI_new[["numgeneR1B"]] = CI_new[["numgeneR1B"]].fillna(0)
  
  weird_genes = ["SNORA","RP11","RP4-","SCARNA","DLEU2","SNORD","CTSLP2"]
  for weird_gene in weird_genes:
    for suff in ["A","B"]:
      ind = CI_new[((CI_new["numgeneR1" + suff] > 2) & (CI_new["geneR1" + suff].str.contains(weird_gene,na=False))) | ((CI_new["numgeneR1" + suff] > 1) & ~(CI_new["gene_strandR1" + suff] == "?") & (CI_new["geneR1" + suff].str.contains(weird_gene,na=False)))].index
      CI_new.loc[ind,"geneR1" + suff] = CI_new.loc[ind,"geneR1" + suff].str.replace("{}[^,]*[,]".format(weird_gene),"",regex=True).str.replace(",{}.*".format(weird_gene),"")
      CI_new.loc[ind,"numgeneR1" + suff] = CI_new.loc[ind,"geneR1" + suff].str.split(",").str.len()#.astype("Int32")
  CI_new["shared_gene"] = [",".join([x for x in a.split(",") if x in b.split(",")]) for a,b in zip(CI_new["geneR1A"],CI_new["geneR1B"])]
  
  CI_new["num_shared_genes"] = CI_new["shared_gene"].str.split(",").str.len()
  CI_new.loc[CI_new["shared_gene"] == "","num_shared_genes"] = 0
  ind = CI_new[(CI_new["num_shared_genes"] > 0) & ((CI_new["numgeneR1A"] > 1) | (CI_new["numgeneR1B"] > 1))].index
  CI_new.loc[ind,"geneR1A"] = CI_new.loc[ind]["shared_gene"].str.split(",").str[-1]
  CI_new.loc[ind,"geneR1B"] = CI_new.loc[ind]["shared_gene"].str.split(",").str[-1]
  CI_new["geneR1A_uniq"] = CI_new["geneR1A"]
  CI_new["geneR1B_uniq"] = CI_new["geneR1B"]
  
  ind = CI_new[(CI_new["numgeneR1A"] > 1) & (CI_new["num_shared_genes"] == 0)].index
  CI_new.loc[ind,"geneR1A_uniq"] = CI_new.loc[ind]["geneR1A"].str.split(",").str[-1]
  ind = CI_new[(CI_new["numgeneR1B"] > 1) & (CI_new["num_shared_genes"] == 0)].index
  CI_new.loc[ind,"geneR1B_uniq"] = CI_new.loc[ind]["geneR1B"].str.split(",").str[-1]
  for let in ["A","B"]:
  
    CI_new = CI_new.merge(gene_strand_info,how="left",left_on = ["geneR1{}_uniq".format(let),"chrR1{}".format(let)], right_on=["gene_name","seqname"])
    CI_new = CI_new.rename(columns={"strand" : "gene_strandR1{}_new".format(let)})
    CI_new = CI_new.drop(["gene_name","seqname"],axis=1)

  # if the library is stranded, we want to keep the read strand; the genes should all come from that strand as well (when not, it seems to be due to strand ambiguity, i.e. the gene appears on both strands)
  if stranded_library:
    for let in ["A","B"]:
      CI_new["gene_strandR1{}_new".format(let)] = CI_new["read_strandR1{}".format(let)]

  ind = CI_new[((((CI_new["gene_strandR1A_new"] != CI_new["read_strandR1A"]) & (CI_new["gene_strandR1B_new"] == CI_new["read_strandR1B"])) | ((CI_new["gene_strandR1A_new"] == CI_new["read_strandR1A"]) & (~CI_new["gene_strandR1B_new"].isna()) & (CI_new["gene_strandR1B_new"] != CI_new["read_strandR1B"]))) & (CI_new["gene_strandR1A"] == "?") & (CI_new["num_shared_genes"] == 0)) & (CI_new["numgeneR1A"] > 1)].index
  CI_new.loc[ind,"geneR1A_uniq"] = CI_new.loc[ind]["geneR1A"].str.split(",").str[-2]
  CI_new = CI_new.drop(["gene_strandR1A_new"],axis=1)
  CI_new = CI_new.merge(gene_strand_info,how="left",left_on = ["geneR1A_uniq","chrR1A"], right_on=["gene_name","seqname"])
  CI_new = CI_new.drop(["gene_name","seqname"],axis=1)
  CI_new = CI_new.rename(columns={"strand" : "gene_strandR1A_new"})

  ind = CI_new[(((CI_new["gene_strandR1A_new"] != CI_new["read_strandR1A"]) & (~CI_new["gene_strandR1A_new"].isna())  & (CI_new["gene_strandR1B_new"] == CI_new["read_strandR1B"])) | ((CI_new["gene_strandR1A_new"] == CI_new["read_strandR1A"]) & (CI_new["gene_strandR1B_new"] != CI_new["read_strandR1B"]))) & (CI_new["gene_strandR1B"] == "?") & (CI_new["num_shared_genes"] == 0) & (CI_new["numgeneR1B"] > 1)].index
  CI_new.loc[ind,"geneR1B_uniq"] = CI_new.loc[ind]["geneR1B"].str.split(",").str[-2]
  CI_new = CI_new.drop(["gene_strandR1B_new"],axis=1)
  CI_new = CI_new.merge(gene_strand_info,how="left",left_on = ["geneR1B_uniq","chrR1B"], right_on=["gene_name","seqname"])
  CI_new = CI_new.rename(columns={"strand" : "gene_strandR1B_new"})
  CI_new = CI_new.drop(["gene_name","seqname"],axis=1)

  if stranded_library:
    for let in ["A","B"]:
      CI_new["gene_strandR1{}_new".format(let)] = CI_new["read_strandR1{}".format(let)]
  
  reverse = {"+" : "-", "-" : "+"}
  same = {"-" : "-", "+" : "+"}
  
  ind = CI_new[(CI_new["gene_strandR1B_new"].isna()) & (CI_new["gene_strandR1A_new"] == CI_new["read_strandR1A"])].index
  CI_new.loc[ind,"gene_strandR1B_new"] = CI_new.loc[ind]["read_strandR1B"].map(same)
  
  ind = CI_new[(CI_new["gene_strandR1B_new"].isna()) & (CI_new["gene_strandR1A_new"] != CI_new["read_strandR1A"]) & (~CI_new["gene_strandR1A_new"].isna())].index
  CI_new.loc[ind,"gene_strandR1B_new"] = CI_new.loc[ind]["read_strandR1B"].map(reverse)
  
  ind = CI_new[(CI_new["gene_strandR1A_new"].isna()) & (CI_new["gene_strandR1B_new"] == CI_new["read_strandR1B"])].index
  CI_new.loc[ind,"gene_strandR1A_new"] = CI_new.loc[ind]["read_strandR1A"].map(same)
  
  ind = CI_new[(CI_new["gene_strandR1A_new"].isna()) & (CI_new["gene_strandR1B_new"] != CI_new["read_strandR1B"]) & (~CI_new["gene_strandR1B_new"].isna())].index
  CI_new.loc[ind,"gene_strandR1A_new"] = CI_new.loc[ind]["read_strandR1A"].map(reverse)
  CI_new["refName_newR1"] = ""
  CI_new["geneR1B_uniq"].fillna("",inplace=True)
  CI_new["geneR1A_uniq"].fillna("",inplace=True)
  
  CI_new["reverse"] = False
  ind = CI_new[(CI_new["fileTypeR1"] == "Aligned") & (CI_new["gene_strandR1A_new"] == "-") & (CI_new["juncPosR1A"] < CI_new["juncPosR1B"])].index
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["chrR1B"] + ":" + CI_new.loc[ind]["geneR1B_uniq"].astype(str) + ":" + CI_new.loc[ind]["juncPosR1B"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1B_new"] + "|" + CI_new.loc[ind]["chrR1A"] + ":" + CI_new.loc[ind]["geneR1A_uniq"].astype(str) + ":" + CI_new.loc[ind]["juncPosR1A"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1A_new"]
  CI_new.loc[ind,"reverse"] = True
  name_swap = {}
  for c in CI_new.columns:
    if "R1A" in c:
      name_swap[c] = c.replace("R1A","R1B")
      name_swap[c.replace("R1A","R1B")] = c
  

  if swap_names:
    CI_new.loc[ind] = CI_new.loc[ind].rename(columns=name_swap)
  ind = CI_new[(CI_new["fileTypeR1"] == "Aligned") & (CI_new["gene_strandR1A_new"] == "+")].index
  
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["chrR1A"] + ":" + CI_new.loc[ind]["geneR1A_uniq"] + ":" + CI_new.loc[ind]["juncPosR1A"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1A_new"] + "|" +  CI_new.loc[ind]["chrR1B"] + ":" + CI_new.loc[ind]["geneR1B_uniq"] + ":" + CI_new.loc[ind]["juncPosR1B"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1B_new"]
  
  
  ind = CI_new[(CI_new["fileTypeR1"] == "Chimeric") & (CI_new["gene_strandR1A_new"] != CI_new["read_strandR1A"]) & (CI_new["gene_strandR1B_new"] != CI_new["read_strandR1B"])].index
  
  
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["chrR1B"] + ":" + CI_new.loc[ind]["geneR1B_uniq"] + ":" + CI_new.loc[ind]["juncPosR1B"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1B_new"] + "|" + CI_new.loc[ind]["chrR1A"] + ":" + CI_new.loc[ind]["geneR1A_uniq"] + ":" + CI_new.loc[ind]["juncPosR1A"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1A_new"]
  CI_new.loc[ind,"reverse"] = True
  if swap_names:
    CI_new.loc[ind] = CI_new.loc[ind].rename(columns=name_swap)
  
  ind = CI_new[(CI_new["fileTypeR1"] == "Chimeric") & ((CI_new["gene_strandR1A_new"] == CI_new["read_strandR1A"]) | (CI_new["gene_strandR1B_new"] == CI_new["read_strandR1B"]))].index
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["chrR1A"] + ":" + CI_new.loc[ind]["geneR1A_uniq"].astype(str) + ":" + CI_new.loc[ind]["juncPosR1A"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1A_new"] + "|" +  CI_new.loc[ind]["chrR1B"] + ":" + CI_new.loc[ind]["geneR1B_uniq"].astype(str) + ":" + CI_new.loc[ind]["juncPosR1B"].astype(str) + ":" + CI_new.loc[ind]["gene_strandR1B_new"]
 
  ind1 = CI_new[(CI_new["refName_newR1"] == "") | (CI_new["refName_newR1"].isna())].index # this ind1 is used to simply replace refName_newR1 with the refName_ABR1 

  CI_new.loc[ind1,"refName_newR1"] = CI_new.loc[ind1]["refName_ABR1"]
  ref_dict = pd.Series(CI_new.refName_newR1.values,index=CI_new.refName_ABR1).to_dict()
  rev_dict = pd.Series(CI_new.reverse.values,index=CI_new.refName_ABR1).to_dict()
  CI["refName_newR1"] = CI["refName_ABR1"].map(ref_dict)
  CI["reverse"] = CI["refName_ABR1"].map(rev_dict)
  ind = CI["reverse"].index
  name_swap = {}

  for c in CI.columns:
    if str(CI[c].dtype)[0] == "i":
      CI[c] = CI[c].astype("float32")


    if "R1A" in c:
      name_swap[c] = c.replace("R1A","R1B")
      name_swap[c.replace("R1A","R1B")] = c


  if swap_names:
    CI.loc[ind] = CI.loc[ind].rename(columns=name_swap)
  CI_new = CI

  CI_new["gene_strandR1A"] = CI_new["refName_newR1"].str.split("|").str[0].str.split(":").str[-1]
  CI_new["gene_strandR1B"] = CI_new["refName_newR1"].str.split("|").str[1].str.split(":").str[-1]
  CI_new["juncPosR1A"] = CI_new["refName_newR1"].str.split("|").str[0].str.split(":").str[2].astype("int")
  CI_new["juncPosR1B"] = CI_new["refName_newR1"].str.split("|").str[1].str.split(":").str[2].astype("int")
  CI_new["chrR1A"] = CI_new["refName_newR1"].str.split("|").str[0].str.split(":").str[0]
  CI_new["chrR1B"] = CI_new["refName_newR1"].str.split("|").str[1].str.split(":").str[0]
  CI_new["geneR1A_uniq"] = CI_new["refName_newR1"].str.split("|").str[0].str.split(":").str[1]
  CI_new["geneR1B_uniq"] = CI_new["refName_newR1"].str.split("|").str[1].str.split(":").str[1]
  CI_new.drop("reverse",axis=1,inplace=True)

## adding the junction type to refName_newR1
  CI_new["junc_type"] = ""
  ind = CI_new[(CI_new["chrR1A"] != CI_new["chrR1B"]) | ((CI_new["read_strandR1A"] == CI_new["read_strandR1B"]) & (abs(CI_new["juncPosR1A"] - CI_new["juncPosR1B"])>=1000000) ) & (CI_new["fileTypeR1"]=="Chimeric")].index
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["refName_newR1"] + "|fus"
  CI_new.loc[ind,"junc_type"] = "fus"

  ind = CI_new[(CI_new["chrR1A"] == CI_new["chrR1B"]) & (CI_new["read_strandR1A"] != CI_new["read_strandR1B"])  & (CI_new["fileTypeR1"]=="Chimeric")].index
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["refName_newR1"] + "|sc"
  CI_new.loc[ind,"junc_type"] = "sc"

  ind = CI_new[(CI_new["junc_type"] != "sc") & (CI_new["junc_type"] != "fus") &  ( ((CI_new["read_strandR1A"] == "+") & (CI_new["juncPosR1A"] > CI_new["juncPosR1B"])) | ((CI_new["read_strandR1A"] == "-") & (CI_new["juncPosR1A"] < CI_new["juncPosR1B"])) ) & (CI_new["fileTypeR1"]=="Chimeric")].index
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["refName_newR1"] + "|rev"
  CI_new.loc[ind,"junc_type"] = "rev"

  ind = CI_new[((CI_new["junc_type"] != "sc") & (CI_new["junc_type"] != "fus") & (CI_new["junc_type"] != "rev")) | (CI_new["fileTypeR1"]=="Aligned")].index
  CI_new.loc[ind,"refName_newR1"] = CI_new.loc[ind]["refName_newR1"] + "|lin"
 
  CI_new = CI_new.drop("junc_type", axis=1)
  
  return CI_new


def get_loc_flag(row):
  if "|fus" in row["refName_ABR1"]:
    return 1
  if "|sc" in row["refName_ABR1"]:
    return 0
  r_strand = row["read_strandR1A"]
  if math.isnan(row["juncPosR2A"]):
    return 0
  if "|lin" in row["refName_ABR1"]:
    if r_strand == "+":
      if row["juncPosR2A"] >= row["juncPosR1A"]:
        return 1
      else:
        return 0
    if r_strand == "-":
      if (row["juncPosR2A"]) <= row["juncPosR1A"]:
        return 1
      else:
        return 0
  if "|rev" in row["refName_ABR1"]:
    if (row["juncPosR1A"] <= row["juncPosR2A"] <= row["juncPosR1B"]) or (row["juncPosR1B"] <= row["juncPosR2A"] <= row["juncPosR1A"]):
      return 1
    else:
      return 0
  return -1

def parse_cigar(cigar):
  matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
  val = 0
  for m in matches:
    if m[1] == "M":
      val += int(m[0])
    elif m[1] == "N":
      val += int(m[0])
    elif m[1] == "D":
      val += int(m[0])
  return val

def read_strand(flag, fill_char = np.nan):
  if flag == fill_char:
    return flag
  sign_dict = {"0" : "+", "1" : "-"}
  return sign_dict['{0:012b}'.format(flag)[7]]

def chim_refName(flags, cigars, offsets, rnames, ann, stranded_library):
    sign_dict = {"0" : "+", "1" : "-"}
    signs = []
    for i in range(len(flags)):
        f = flags[i]
        signs.append(sign_dict['{0:012b}'.format(f)[7]])


    # determined by comparing Chimeric.out.junction and Chimeric.out.sam
    if signs[0] == "+":
      cig_val = parse_cigar(cigars[0])
      posFirst = int(offsets[0]) + cig_val - 1
    elif signs[0] == "-":
      posFirst = int(offsets[0])
    else:
      print("flags", flags)

    if signs[1] == "+":
      posSecond = int(offsets[1])
    elif signs[1] == "-":
      cig_val = parse_cigar(cigars[1])
      posSecond = int(offsets[1]) + cig_val - 1



    gene1, strand1 =  ann.get_name_given_locus(rnames[0], posFirst, signs[0], stranded_library)
    gene2, strand2 =  ann.get_name_given_locus(rnames[1], posSecond, signs[1], stranded_library)

    if rnames[0] != rnames[1]:
        juncType = "fus"
    elif signs[0] != signs[1]:
      juncType = "sc"
    elif (signs[0] == "+" and posFirst > posSecond) or (signs[0] == "-" and posFirst < posSecond):
        juncType = "rev"
    elif (signs[0] == "+" and posFirst < posSecond) or (signs[0] == "-" and posFirst > posSecond):
         juncType = "lin"
    else:
        juncType = "err"
    unchanged = "{}:{}:{}:{}|{}:{}:{}:{}|{}".format(rnames[0], gene1, posFirst, strand1, rnames[1], gene2, posSecond, strand2, juncType)

    return unchanged, rnames[0], gene1, int(posFirst), rnames[1], gene2, int(posSecond)


def readObj_refname(read_strand, cigar, seqname, position, ann, fill_char, stranded_library):
  #flag_dict = {0 : "+", 256 : "+", 16 : "-", 272 : "-"}
  #read_strand = flag_dict[flag]

  if "N" not in cigar:
    gene, strand = ann.get_name_given_locus(seqname, int(position), read_strand, stranded_library)
    return "{}:{}:{}".format(seqname,gene,strand), seqname,gene, position, fill_char, fill_char, fill_char

  matches = re.findall(r'(\d+)([A-Z]{1})', cigar)

  # find the largest N (the one to split on)
  max_N_ind = None
  max_N_val = 0
  for i in range(len(matches)):
      m = matches[i]
      if m[1] == "N":
          if int(m[0]) > max_N_val:
              max_N_ind = i
              max_N_val = int(m[0])

  # get the first base of the junction
  offset1 = position
  for i in range(max_N_ind):
      m = matches[i]
      if m[1] in ["M","N","D"]:
          offset1 += int(m[0])

  # get the second base of the junction
  offset2 = offset1 + max_N_val
  for i in range(max_N_ind + 1, len(matches) + 1):
      m = matches[i]
      if m[1] == "M":
          break

      elif m[1] in ["N","D"]:
          offset2 += int(m[0])
  offset1 -= 1
  gene1, strand1 = ann.get_name_given_locus(seqname, offset1, read_strand, stranded_library)
  gene2, strand2 = ann.get_name_given_locus(seqname, offset2, read_strand, stranded_library)

  read_class = "lin"

  return "{}:{}:{}:{}|{}:{}:{}:{}|{}".format(seqname, gene1, int(offset1), strand1,seqname, gene2, int(offset2), strand2, read_class), seqname, gene1, offset1, seqname, gene2, offset2


def get_SM(cigar, fill_char = np.nan):
  if not isinstance(cigar, str):
    return fill_char, fill_char, fill_char, fill_char
  matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
  M = 0
  S = 0
  I = 0
  D = 0
  for m in matches:
    if m[1] == "M":
      M += int(m[0])
    elif m[1] == "S":
      S += int(m[0])
    elif m[1] == "I":
      I += int(m[0])
    elif m[1] == "D":
      D += int(m[0])

  return M, S, I, D

def entropy(kmer, c=5):
    """Calculate the entropy of a kmer using cmers as the pieces of information"""
    num_cmers = len(kmer) - c + 1
    cmers = []
    for i in range(num_cmers):
        cmers.append(kmer[i:i+c])
    Ent = 0
    for cmer in set(cmers):

        prob = cmers.count(cmer)/float(num_cmers)
        Ent -= prob*np.log(prob)
    return Ent

def count_stretch(s):
  counts = {"A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0}
  for i in range(len(s)):
    if i == 0:
      curr_stretch = 1
    elif s[i] != s[i-1]:
      counts[s[i-1]] = max(counts[s[i-1]], curr_stretch)
      curr_stretch = 1
    else:
      curr_stretch += 1
  counts[s[-1]] = max(counts[s[-1]], curr_stretch)
  return counts

def nmm(MD):
  return len(''.join(filter(["A","C","G","T"].__contains__, MD)))

def split_cigar_align(cigar, fill_char = np.nan):
  if "N" not in cigar:
    return cigar, fill_char
  matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
  
  # find the largest N (the one to split on)
  max_N_ind = None
  max_N_val = 0
  for i in range(len(matches)):
      m = matches[i]
      if m[1] == "N":
          if int(m[0]) > max_N_val:
              max_N_ind = i
              max_N_val = int(m[0])
  cigar1 = "".join(["".join(m) for m in matches[:max_N_ind]])
  cigar2 = "".join(["".join(m) for m in matches[max_N_ind + 1:]])
  return cigar1, cigar2

def split_cigar_chim(cigar):
    matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
    if (matches[0][1] in ["S","H"]) and (not(matches[-1][1] in ["S","H"])):
      cigar1 = "".join(["".join(m) for m in matches[1:]])
    elif (matches[-1][1] in ["S","H"]) and (not(matches[0][1] in ["S","H"])):
      cigar1 = "".join(["".join(m) for m in matches[:-1]])
    else:
      assert (matches[0][1] in ["S","H"]) and (matches[-1][1] in ["S","H"])
      if int(matches[0][0]) > int(matches[-1][0]):
        cigar1 = "".join(["".join(m) for m in matches[1:]])
      else:
        cigar1 = "".join(["".join(m) for m in matches[:-1]])
    return cigar1