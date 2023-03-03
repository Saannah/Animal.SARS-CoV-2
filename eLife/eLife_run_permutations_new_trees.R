#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ape)
library(stringr)
library(phylotools)
library(ggtree)
library(tidyverse)
library(ggtree)
library(tools)
library(mvSLOUCH) #for fitch.mvsl()

find_intro = function(species, number){
  fromto = c("human", "animal")
  print(number)
  #tree = read.tree("deer_subtree_5_trimmed.nwk")
  #tree = read.tree(paste0(species,"_subtree_",number,"_trimmed.nwk"))
  ## shuffle tips
  tree$tip.label = tree$tip.label[sample(1 : length(tree$tip.label))]
  
  
  
  trait = rep("human", length(tree$tip.label))
  trait[which(str_detect(tree$tip.label, paste0("/", species, "/")))] = "animal"
  labels = data.frame(id = tree$tip.label, trait = trait)
  
  
  method = "ml"
  ml_asr_model = "ER"
  root_state = NULL
  transition_fromto = fromto
  
  
  ################################# data loading and cleaning starts here
  strainlist = tree$tip.label #list of all strain id's we have
  labels= labels[which(labels[,1] %in% strainlist),]
  ################################# data loading and cleaning ends here
  
  
  ###### check if tree is rooted
  if (!is.rooted(tree)){
    tree = root(tree, outgroup = 1, resolve.root = TRUE)
  }
  
  ###### add an attribute to the tree with its bootstrap support value 
  tree$bootstrap = "-1"
  for(i in 2:length(tree$node.label)){
    tree$bootstrap[i] = unlist(str_split(tree$node.label[i], "/"))[2]
  }
  tree$bootstrap = as.numeric(tree$bootstrap)
  ################################# converting the tree to a fully dichotomous tree, starts here
  dichotomous_tree = multi2di(tree, random = FALSE)   #format change of tree
  dichotomous_tree$edge.length[dichotomous_tree$edge.length<=0] = 1e-8  #set very small values to zero and negative values (needed for ace to run properly)
  ################################# converting the tree to a fully dichotomous tree, starts here
  
  
  
  ################################# creating the labels for reconstruction, starts here
  # "labels" is the matrix that includes the tip labels for the tree
  # the dimension of labels is ntips x 2. the first column has sample IDs and the second column has the states corresponding to each sample.
  
  states = unique(labels[,2]) #the unique values of the states DO NOT CHANGE THE ORDER OF THIS VECTOR IN THE REMAINDER OF THE CODE
  n_states = length(states)
  
  
  mpr_tipstates = rep(0, length(dichotomous_tree$tip.label))
  
  for (i in 1:n_states){ #loop to assign the tip labels for the traits 
    s = states[i]
    whos = labels[which(labels[,2] == s),1]
    mpr_tipstates[match(whos, tree$tip.label)] = i
  }
  
  dichotomous_tree$tip.label -> names(mpr_tipstates)
  ################################# creating the labels for reconstruction, ends here
  
  
  
  ################################# MPR reconstruction starts here
  # Whole-tree MPR reconstruction, fitch.mvsl function is used from "mvslouch" package.
  # both 'deltran' and 'acctran' settinge will be applied.
  
  if (method == "acc"){
    
    if(is.null(root_state)){mpr_recon_acctran = fitch.mvsl(dichotomous_tree, mpr_tipstates, deltran = FALSE, acctran = TRUE)}
    if(!is.null(root_state)){mpr_recon_acctran = fitch.mvsl(dichotomous_tree, mpr_tipstates, deltran = FALSE, acctran = TRUE, root = root_state)}
  }
  
  if (method == "del"){
    
    if(is.null(root_state)){mpr_recon_deltran = fitch.mvsl(dichotomous_tree, mpr_tipstates, deltran = TRUE, acctran = FALSE)}
    if(!is.null(root_state)){mpr_recon_deltran = fitch.mvsl(dichotomous_tree, mpr_tipstates, deltran = TRUE, acctran = FALSE, root = root_state)}
  }
  
  
  ################################# MPR reconstruction ends here
  
  
  ################################# ML reconstruction ends here
  # labels vector format is same as the MPR analysis, the same vector was used
  
  if (method == "ml"){
    ml_recon = ace(mpr_tipstates, dichotomous_tree, type = "discrete", model="ER")
  }
  #
  ################################# ML reconstruction ends here
  
  
  #### NOTE:
  # the output of the mpr function corresponds to the edge attribute of thhe tree
  # the output of the ace fucntion increments from the top of the tree (starting from the root)
  
  
  
  comparison = comparePhylo(tree, dichotomous_tree) #compares the original polytomous tree to the dichotomous tree
  if(is.null(comparison$NODES)){
    e1 = length(tree$tip.label)+1
    e2 = (length(tree$tip.label) + (tree$Nnode))
    ns = seq(e1, e2)
    seq = paste0("(", as.character(ns), ")")
    comparison$NODES = data.frame(tree = seq, dichotomous_tree = seq)
  }
  
  comp = matrix(0, dim(comparison$NODES)[1], dim(comparison$NODES)[2])
  
  
  for (i in 1:dim(comparison$NODES)[1]){ #loop for extracing the node number out of the compare$node table returned by comparePhylo()
    for (j in 1:dim(comparison$NODES)[2]){
      string = comparison$NODES[i,j]
      str = stringr::str_extract(string = comparison$NODES[i,j], pattern = "(?<=\\().*(?=\\))")
      comp[i,j] = as.numeric(str) #replacing the string number to numeric value, optional, might have to change later
    }
  }
  
  
  
  ################################# cleaning up ace output starts here
  
  # the columns of ml_recon$lik.anc have the same name as tip state numbers. meaning that column names will be numbers from 1 to "n_state"
  
  
  if (method == "ml"){
    
    
    ml_asr = ml_recon$lik.anc #the matrix containing likelihoods for everynode
    ml_colnames = colnames(ml_asr)
    
    ace_dataframe = matrix(0, dim(comp)[1], 2)
    lik_animal = vector()
    for (i in 1:dim(comp)[1]){
      row = comp[i,2] - length(dichotomous_tree$tip.label)
      s = which.max(ml_asr[row,]) #the states are numbered from 1 to n, so the column number is the state itself
      ace_dataframe[i,1] = i + length(tree$tip.label) #node number in the actual tree
      ace_dataframe[i,2] = s #the state of that node
      lik_animal = c(lik_animal, ml_asr[row,2])
    }
    
    ace_dataframe = data.frame(node = ace_dataframe[,1], state = ace_dataframe[,2])
  }
  ################################# cleaning up ace output ends here
  
  
  
  
  ################################# cleaning up MPR output starts here
  
  
  
  if (method == "del"){
    ## Deltran
    mpr_polytomous = matrix(0, dim(comp)[1], 2)
    mpr_polytomous[1,2] = as.numeric(mpr_recon_deltran$root_regime) #assigning root state
    mpr_polytomous[1,1] = as.numeric(comp[1,2])
    for (i in 2:dim(comp)[1]){
      n_di = comp[i,2]
      in_edge = which(dichotomous_tree$edge[,2] == n_di)
      mpr_polytomous[i,2] = as.numeric(mpr_recon_deltran$branch_regimes[in_edge])
      mpr_polytomous[i,1] = as.numeric(comp[i,1])
    }
    del_dataframe = data.frame(node = mpr_polytomous[,1], state = mpr_polytomous[,2])
  }
  
  if (method == "acc"){
    ## acctran
    mpr_polytomous = matrix(0, dim(comp)[1], 2)
    mpr_polytomous[1,2] = as.numeric(mpr_recon_acctran$root_regime) #assigning root state
    mpr_polytomous[1,1] = as.numeric(comp[1,2])
    for (i in 2:dim(comp)[1]){
      n_di = comp[i,2]
      in_edge = which(dichotomous_tree$edge[,2] == n_di)
      mpr_polytomous[i,2] = as.numeric(mpr_recon_acctran$branch_regimes[in_edge])
      mpr_polytomous[i,1] = as.numeric(comp[i,1])
    }
    acc_dataframe = data.frame(node = mpr_polytomous[,1], state = mpr_polytomous[,2])
  }
  
  
  # converting transition_fromto to numbers [1 to n_states]
  # THE ORDER OF THE "states" VARIABLE MATTERS!!!
  from = transition_fromto[1]
  to = transition_fromto[2]
  n_from = which(states == from) #numerical value for the given states (from)
  n_to = which(states == to) #numerical value for the given states (to)
  
  
  
  ########## transitions:
  
  
  ################################# transition nodes starts here
  
  #i put this conditions to later be able to manipulate data types that are different in type among different methods, might remove later
  if (method == "ml"){
    all = data.frame(node = c(1:length(tree$tip.label), ace_dataframe$node), state = c(mpr_tipstates, ace_dataframe$state))
  }
  
  if (method == "acc"){
    all = data.frame(node = c(1:length(tree$tip.label), acc_dataframe$node), state = c(mpr_tipstates, acc_dataframe$state))
  }
  
  if (method == "del"){
    all = data.frame(node = c(1:length(tree$tip.label), del_dataframe$node), state = c(mpr_tipstates, del_dataframe$state))
  }
  
  ## add a column to all dataframe
  all$labels = states[1]
  all$labels[which(all$state == 2)] = states[2]
  
  #sort edge table by the second column (children)
  edge = tree$edge
  edge = data.frame(parent = edge[,1], child = edge[,2])
  edge = edge[order(edge$child),]
  
  
  # finding *all* transition events in the tree
  o = vector()
  for (i in 1:dim(edge)[1]){
    s = all[edge$child[i], 2]
    if(s == n_to){
      o = append(o,i)
    }
  }
  transitions_childeren = edge[o,]
  z = vector()
  for (i in 1:dim(transitions_childeren)[1]){
    s = all[transitions_childeren$parent[i], 2]
    if(s == n_from){
      z = append(z,i)
    }
  }
  transitions_parents = transitions_childeren[z,]
  #trans_nodes includes ALL events where e non-QC node has a QC child. this includes embedded transitions and singletons 
  trans_nodes = sort(unique(transitions_parents$child)) 
  trans_nodes_store = trans_nodes
  
  
  
  
  if(length(trans_nodes) == 0){
    out_obj = list(unfiltered_MRCAin = 0,
                   unfiltered_MRCAout = 0)
    return(out_obj)
  }
  
  
  #function for finding the parents of a node up to the root
  # inputs: tree, node number outputs: a vector of all parents of the node up to the root
  parent_finder = function(tree, nodes){ #function starts here
    parents_list = list()
    for (node in nodes){
      edge = tree$edge
      parents = vector()
      flag = 0
      root = length(tree$tip.label) + 1
      
      parents = append(parents, node)
      while(flag == 0){
        if (node == root){
          flag = 1
          next
        }
        ind = which(edge[,2] == node)
        node = edge[ind, 1]
        parents = append(parents, node)
        if (node == root){
          flag = 1
        }
      }
      parents_list = append(parents_list, list(parents))
    }
    return(parents_list)
  } #function ends here
  
  parents_list = parent_finder(tree, trans_nodes)
  
  #function for finding transition parents of transition nodes
  # this function gets a list of vectors representing nodes ancestors, returns transision events among them
  
  
  checked_list = vector()
  transitions_obj = list()
  for (i in 1:length(parents_list)){ #loop starts here
    p = unlist(parents_list[i])
    if (p[1] %in% checked_list){
      next
    }
    
    parents = vector()
    for (j in 1:length(p)){
      if (p[j] %in% trans_nodes){
        parents = append(parents, p[j])
      }
    }
    checked_list = append(checked_list, parents)
    transitions_obj = append(transitions_obj, list(parents))
    
  }#loop ends here
  
  independent_transitions = vector()
  for (i in 1:length(transitions_obj)){ #loop starts here
    p = unlist(transitions_obj[i])
    independent_transitions = append(independent_transitions, p[length(p)])
  }#loop ends here
  independent_transitions = unique(independent_transitions)
  ################################# transition nodes ends here
  
  
  new_independent_transitions = independent_transitions
  new_trans_nodes = trans_nodes
  
  
  
  
  # finding the non-qc parents of "new transitions"
  new_independent_transitions_parents = vector()
  for (i in 1:length(new_independent_transitions)){
    ind = which(edge[,2] == new_independent_transitions[i])
    new_independent_transitions_parents = append(new_independent_transitions_parents, edge[ind,1])
  }
  
  MRCAout = list(new_independent_transitions_parents)
  MRCAin = list(new_independent_transitions)
  
  intros = list(MRCAout = MRCAout, MRCAin = MRCAin, reconstruction = all)
  
  
  
  # plot_intros(MRCAin)
  # print("mrca")
  # print(MRCAin)
  out_obj = length(unlist(MRCAin))
  return(out_obj)}


species = args[1]
subtree = args[2]
print("species:")
print(species)
print("subtree:")
print(subtree)


#load tree
path = paste0(species, "_subtree_", subtree, "_trimmed.nwk")
tree = read.tree(path)


#tic()
tab = mapply(find_intro, species, 1:10)
#toc()

tab = as.vector(t(tab))
#write file
write.table(tab, file = paste0(species, "_", subtree, "permutations_human2animal.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)