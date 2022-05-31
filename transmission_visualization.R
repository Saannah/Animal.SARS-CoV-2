library(ape)
library(stringr)
library(phylotools)
library(ggtree)


tree = read.tree("/Users/sana/Documents/Jan2022/plots_18032022/animal_to_human/cat_subtree_1_trimmed.nwk")
species = "cat"
number = 1
fromto = c("animal", "human")

trait = rep("human", length(tree$tip.label))
trait[which(str_detect(tree$tip.label, "/cat/"))] = "animal"
labels = data.frame(id = tree$tip.label, trait = trait)

find_intro = function(tree, metadata, species, number, fromto){
  ## intros for animal sars
  ### load libraries
  #setwd("/Users/sana/Documents/Jan2022/split_trees/mink/")
  
  
  method = "ml"
  ml_asr_model = "ER"
  root_state = NULL
  #fromto = c("animal","human")
  
  
  
  
  
  ### filter metadata based on tips existing in tree:
  #metadata = metadata[which(metadata[, 23] %in% tree$tip.label),]
  #trait = rep("animal", length(metadata[, 7]))
  #trait[which(metadata[, 8] == "Human")] = "human"
  
  #labels = data.frame(id = metadata[, 23], trait = trait)
  
  #### specify transition from/to
  transition_fromto = fromto
  
  
  ################################# data loading and cleaning starts here
  strainlist = tree$tip.label #list of all strain id's we have
  labels= labels[which(labels[,1] %in% strainlist),]
  ################################# data loading and cleaning ends here
  
  
  ###### check if tree is rooted
  if (!is.rooted(tree)){
    tree = root(tree, outgroup = 1, resolve.root = TRUE)
  }
  
  
  
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
  
  
  
  ### define plotting function
  plot_intros = function(MRCAin){ # plot intros function starts here
    #MRCAin = vector()
    library(tidyverse)
    library(ggtree)
    library(tools)
    library(mvSLOUCH) #for fitch.mvsl()
    library(ape)
    
    #load metadata
    
    # slashes = str_locate_all(tree$tip.label, "_")
    # tree_tip_lab = vector()
    # for(i in 1:length(slashes)){
    #   if(metadata[which(metadata[,23] == tree$tip.label[i]), 8] == "Human"){
    #     tree_tip_lab[i] = paste0("/","human",substr(tree$tip.label[i], unlist(slashes[i])[1], unlist(slashes[i])[2]))
    #     next
    #   }
    #   tree_tip_lab[i] = substr(tree$tip.label[i], unlist(slashes[i])[1], unlist(slashes[i])[3])
    # }
    # tree_tip_lab = str_replace_all(tree_tip_lab,"_","/")
    # 
    # p_tree = tree
    # p_tree$tip.label = tree_tip_lab
    # p_tree$edge.length[which(p_tree$edge.length <= 0)] = 1e-5
    
    #  #path = paste0("/home/sanna/projects/def-shapiro/sanna/animal_sars/13-subtree_intros/besthits_names/",species,"_besthits.txt")
    #     besthit_names = read.delim("deer_besthits.txt", header = FALSE)
    #     besthit_names = besthit_names[,1]
    #     besthit_names = str_remove_all(besthit_names,">")
    #     besthit_names = str_replace_all(besthit_names,"/","_")
    #     besthit_names = str_replace_all(besthit_names,"\\|","_")
    #     
    #     
    #     trait = rep("animal", length(tree$tip.label))
    #     trait[which(str_detect(p_tree$tip.label,"human"))] = "human"
    #     trait[which(tree$tip.label %in% besthit_names)] = "besthit"
    
    # labels = data.frame(id = metadata[, 23], traits = metadata[, 24])
    # row.names(labels) = tree$tip.label
    # animal_ids = which(labels[,2] == species)
    # besthit_ids = which(labels[,2] == paste0("besthit_", species))
    # labels[,2] = rep("human", dim(labels)[1])
    # labels[animal_ids,2] = "animal"
    # labels[besthit_ids,2] = "besthit"
    # 
    # #### CREATING PIES DF FOR PLOTTING
    # #human is state 2, animal is state 1
    # ml_asr = ml_recon$lik.anc #the matrix containing likelihoods for everynode
    # ml_colnames = colnames(ml_asr)
    # 
    # ace_dataframe = matrix(0, dim(comp)[1], 3)
    # 
    # for (i in 1:dim(comp)[1]){
    #   row = comp[i,2] - length(dichotomous_tree$tip.label)
    #   s = ml_asr[row,] #the states are numbered from 1 to n, so the column number is the state itself
    #   ace_dataframe[i,1] = i + length(tree$tip.label) #node number in the actual tree
    #   ace_dataframe[i,2:3] = s #the state of that node
    # }
    # 
    # ace_dataframe = data.frame(node = ace_dataframe[,1], animal = ace_dataframe[,2], human = ace_dataframe[,3])
    
    #pies = nodepie(ace_dataframe, cols=2:3, color=c("orange", "cadetblue3"), alpha=0.8)
    #node = [(length(tree$tip.label)+1):(length(tree$tip.label)+length(tree$))]
    
    #1 for human 2 for animal (check this from the order of states)
    df = ace_dataframe
    node_label = rep(states[1], length(df$node))
    node_label[which(df$state == 2)] = states[2]
    node_lab = data.frame(node = df$node, trait = node_label)
    
    
    #  node_df = all[which(all$node > (length(tree$tip.label))),]
    #     node_df$state = tree$node.label
    #     likelihood_df = data.frame(node = ace_dataframe$node, lik = lik_animal)
    #     halves = likelihood_df[which(likelihood_df$lik == 0.5),]
    
    human_nodes = df$node[which(node_label == "human")]
    animal_nodes = df$node[which(node_label == "animal")]
    
    
    p = ggtree(tree) %<+% labels +  scale_color_manual(values=c("orange", "cadetblue3"),labels = c("animal","human")) +
      geom_tippoint(aes(color = trait), size = 0.5) + geom_point2(aes(subset=(node %in% human_nodes)),color="cadetblue3",size=2.5, shape = 18) +
      geom_point2(aes(subset=(node %in% animal_nodes)),color="orange",size=2.5, shape = 18) +
      geom_point2(aes(subset=(node %in% unlist(MRCAin))),color="red",size=2.5, shape = 25) +
      geom_nodelab(size = 0.5, hjust=1.5, vjust=-1.5)
    
    
    
    #geom_text2(aes(subset=(node %in% node_df$node)), label = node_df$state, size=0.5, color = 'blue')
    #geom_text2(aes(subset=(node %in% halves$node)), label = halves$lik, size=2.5, color = 'blue')
    #path = paste0("/home/sanna/projects/def-shapiro/sanna/animal_sars/13-subtree_intros/",species,"_subtree",number,"_intros.pdf")
    ggsave(file =paste0(species, "_", number, ".pdf"), plot = p, height = 8, width  = 10, limitsize = FALSE)
    
  } # plot intros function ends here
  
  
  
  
  #returns an empty vector if there are no intros
  if (length(trans_nodes) == 0){
    v = vector()
    print("no intros")
    plot_intros(v)
    return(0)
    #print("no intros in this direction")
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
  
  
  
  
  
  plot_intros(MRCAin)
  print("mrca")
  print(MRCAin)
  
  if (length(trans_nodes)==0){return(0)}
  return(length(unlist(MRCAin)))
} #function closes here




