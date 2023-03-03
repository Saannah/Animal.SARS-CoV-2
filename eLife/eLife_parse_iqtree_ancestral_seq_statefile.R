#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

filename = args[1]

## load required packages
library(phylotools)
library(stringr)

## load the table
tab=read.table(paste0(filename, ".state"),header=TRUE)
print("table loaded")
# the length of the alignment for this study is: 30580
nodes = unique(tab$Node)

paste_seq = function(node_name){
  sub_tab = tab[which(tab$Node == node_name), ]
  seq_text = paste0(sub_tab$State, collapse = "")
  return(seq_text)
}

msa = lapply(nodes, paste_seq)

#write msa to a fasta
library(seqinr)
write.fasta(names = as.list(nodes), sequences = msa, file.out = paste0(filename, ".ancestral_msa.fasta"))
print("fasta written")

#sanity check: load tree and check if number of nodes is correct:
tree = read.tree(paste0(filename, ".treefile"))
if(length(nodes) == tree$Nnode){
  print("clean exit")
}
