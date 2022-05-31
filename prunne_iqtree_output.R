library(ape)
library(ggtree)
library(stringr)



metadata = read.delim("/project/6006375/sanna/ALL_FILES/post_align/post_align_metadata.tsv")


trim_tree = function(species, num){
print(species)
#load tree
path = paste0("/home/sanna/projects/def-shapiro/sanna/7-subtrees/run_iqtree/", species, "/", species, "_subtree_msa_", num, ".fasta.treefile")
tree = read.tree(path)
print("tree loaded")

tree$tip.label = str_replace_all(tree$tip.label, "___","\\|")
tree$tip.label = str_replace_all(tree$tip.label, "__","/")
print("tree names fixex")

#sanity check
print("sanity check: tree names are ok:")
if(length(which(tree$tip.label %in% metadata[,23])) == length(tree$tip.label)){print("sanity check ok")} else print("names wronng!")

## trim the outlier branches ## option 1: plot histogram and do it approximately
gtree = ggtree(tree)
tab = gtree$data
#histogram = hist(tab$x, nclass=1000, xlim = c(0,0.005))
#trimmed_tree = drop.tip(tree, tab$label[which(tab$x > 0.002)], trim.internal = TRUE)

if(length(boxplot.stats(tab$x)$out) > 0){
threshold = min(boxplot.stats(tab$x)$out[which(boxplot.stats(tab$x)$out > mean(tab$x))])
trimmed_tree = drop.tip(tree, tab$label[which(tab$x >= threshold)], trim.internal = TRUE)
print("tree trimmed.")} else {print("no outliers"); trimmed_tree = tree}

print("length tree:")
print(length(trimmed_tree$tip.label))

#write mew tree
path = paste0("/home/sanna/projects/def-shapiro/sanna/7-subtrees/trim_iqtree_output/", species, "_subtree_", num, "_trimmed.nwk")
write.tree(trimmed_tree, file = path)
print("tree written")
}

species = c("cat", "dog", "mink", "deer")
for(s in species){
   for(i in 1:10){
      print(i)
      trim_tree(s, i)
   }

}