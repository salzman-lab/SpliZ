#!/usr/bin/env Rscript

# Written By : Roozbeh Dehghannasiri (rdehghan@stanford.edu)
# this script takes the permutation file and finds the most variable splize sites for genes with permutation p-value <0.05 
# it finds up to 3 splice sites for each eigenvector (1st, 2nd, 3rd)
#it creates three output files corresponding to the splice sites for each eigenvector 

library(data.table)
library(Rfast)

args <- commandArgs(trailingOnly = TRUE)
p_value_file = args[1]
first_evec_file = args[2]
second_evec_file = args[3]
third_evec_file = args[4]
libraryType = args[5]

p_value = fread(p_value_file,sep="\t",header=TRUE)

## I want to select the top 20 and top 50 genes with FDR < 0.05
if (libraryType == "SS2") {
  p_value = p_value[perm_pval_adj_svd_z0<0.05]
}



topgenes = unique(p_value$gene)
print(paste("number of genes to run",length(topgenes)))

gene_to_plot = c() # I get these vectors to build a data table so that their dot plots can be made automatically
coordinate_to_plot = c()
let_to_plot = c()
for (counter in 1:length(topgenes)){
  gene = topgenes[counter]  # name of the gene
  tryCatch({
    geneMat_file = paste(gene, ".geneMat", sep="")
    loadings = fread(geneMat_file)
    loadings_sq = loadings[1,]^2
    top_site = names(loadings_sq)[loadings_sq==max(loadings_sq)]
    coordinate_to_plot = c(coordinate_to_plot,strsplit(top_site,split = "_")[[1]][1])
    let_to_plot = c(let_to_plot,strsplit(top_site,split = "_")[[1]][2])
    gene_to_plot = c(gene_to_plot,gene)
    
    # I copy for the second and third only if they have at least 10% of loadings
    if (Rfast::nth(as.matrix(loadings_sq), 2, descending = T) > 0.1){
      second_top_site = names(loadings_sq)[loadings_sq == Rfast::nth(as.matrix(loadings_sq), 2, descending = T)]
      coordinate_to_plot = c(coordinate_to_plot,strsplit(second_top_site,split = "_")[[1]][1])
      let_to_plot = c(let_to_plot,strsplit(second_top_site,split = "_")[[1]][2])
      gene_to_plot = c(gene_to_plot,gene)
    }
    if (Rfast::nth(as.matrix(loadings_sq), 3, descending = T) > 0.1){
      third_top_site = names(loadings_sq)[loadings_sq == Rfast::nth(as.matrix(loadings_sq), 3, descending = T)]
      coordinate_to_plot = c(coordinate_to_plot,strsplit(third_top_site,split = "_")[[1]][1])
      let_to_plot = c(let_to_plot,strsplit(third_top_site,split = "_")[[1]][2])
      gene_to_plot = c(gene_to_plot,gene)
    }
    
    top_site = ""
    second_top_site = ""
    third_top_site = ""
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
to_plot = data.table(gene_to_plot,let_to_plot,coordinate_to_plot)
names(to_plot) = c("gene","let","end")

write.table(to_plot, first_evec_file, sep = "\t", row.names = FALSE, quote = FALSE)

##############################
#### second eigen vector #####
##############################

gene_to_plot = c() # I get these vectors to build a data table so that their dot plots can be made automatically
coordinate_to_plot = c()
let_to_plot = c()
for (counter in 1:length(topgenes)){
  gene = topgenes[counter]  # name of the gene
  tryCatch({
    geneMat_file = paste(gene, ".geneMat", sep="")
    loadings = fread(geneMat_file)
    loadings_sq = loadings[2,]^2
    top_site = names(loadings_sq)[loadings_sq==max(loadings_sq)]
    coordinate_to_plot = c(coordinate_to_plot,strsplit(top_site,split = "_")[[1]][1])
    let_to_plot = c(let_to_plot,strsplit(top_site,split = "_")[[1]][2])
    gene_to_plot = c(gene_to_plot,gene)
    
    # I copy for the second and third only if they have at least 10% of loadings
    if (Rfast::nth(as.matrix(loadings_sq), 2, descending = T) > 0.1){
      second_top_site = names(loadings_sq)[loadings_sq == Rfast::nth(as.matrix(loadings_sq), 2, descending = T)]
      coordinate_to_plot = c(coordinate_to_plot,strsplit(second_top_site,split = "_")[[1]][1])
      let_to_plot = c(let_to_plot,strsplit(second_top_site,split = "_")[[1]][2])
      gene_to_plot = c(gene_to_plot,gene)
    }
    if (Rfast::nth(as.matrix(loadings_sq), 3, descending = T) > 0.1){
      third_top_site = names(loadings_sq)[loadings_sq == Rfast::nth(as.matrix(loadings_sq), 3, descending = T)]
      coordinate_to_plot = c(coordinate_to_plot,strsplit(third_top_site,split = "_")[[1]][1])
      let_to_plot = c(let_to_plot,strsplit(third_top_site,split = "_")[[1]][2])
      gene_to_plot = c(gene_to_plot,gene)
    }
    
    top_site = ""
    second_top_site = ""
    third_top_site = ""
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
to_plot = data.table(gene_to_plot,let_to_plot,coordinate_to_plot)
names(to_plot) = c("gene","let","end")

write.table(to_plot, second_evec_file, sep = "\t", row.names = FALSE, quote = FALSE)


##############################
#### third eigen vector #####
##############################

gene_to_plot = c() # I get these vectors to build a data table so that their dot plots can be made automatically
coordinate_to_plot = c()
let_to_plot = c()
for (counter in 1:length(topgenes)){
  gene = topgenes[counter]  # name of the gene
  tryCatch({
    geneMat_file = paste(gene, ".geneMat", sep="")
    loadings = fread(geneMat_file)
    loadings_sq = loadings[3,]^2
    top_site = names(loadings_sq)[loadings_sq==max(loadings_sq)]
    coordinate_to_plot = c(coordinate_to_plot,strsplit(top_site,split = "_")[[1]][1])
    let_to_plot = c(let_to_plot,strsplit(top_site,split = "_")[[1]][2])
    gene_to_plot = c(gene_to_plot,gene)
    
    # I copy for the second and third only if they have at least 10% of loadings
    if (Rfast::nth(as.matrix(loadings_sq), 2, descending = T) > 0.1){
      second_top_site = names(loadings_sq)[loadings_sq == Rfast::nth(as.matrix(loadings_sq), 2, descending = T)]
      coordinate_to_plot = c(coordinate_to_plot,strsplit(second_top_site,split = "_")[[1]][1])
      let_to_plot = c(let_to_plot,strsplit(second_top_site,split = "_")[[1]][2])
      gene_to_plot = c(gene_to_plot,gene)
    }
    if (Rfast::nth(as.matrix(loadings_sq), 3, descending = T) > 0.1){
      third_top_site = names(loadings_sq)[loadings_sq == Rfast::nth(as.matrix(loadings_sq), 3, descending = T)]
      coordinate_to_plot = c(coordinate_to_plot,strsplit(third_top_site,split = "_")[[1]][1])
      let_to_plot = c(let_to_plot,strsplit(third_top_site,split = "_")[[1]][2])
      gene_to_plot = c(gene_to_plot,gene)
    }
    
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
to_plot = data.table(gene_to_plot,let_to_plot,coordinate_to_plot)
names(to_plot) = c("gene","let","end")

write.table(to_plot, third_evec_file, sep = "\t", row.names = FALSE, quote = FALSE)
