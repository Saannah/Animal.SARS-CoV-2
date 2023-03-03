library(ape)
library(stringr)
library(phylotools)
library(ggtree)
library(tidyverse)
library(ggtree)
library(tools)
library(mvSLOUCH) #for fitch.mvsl()

##load and run find_intros_func

################ plot main figure trees ##################
#cat plot subtree 2
cat = find_intro("cat", 2, c("animal", "human"))
cat_plot = ggtree(cat$tree, aes(color = labels, size = labels), layout = 'circular', branch.length = 'none') %<+% cat$all +  
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  scale_size_manual(values = c(0.9, 0.5), labels = c("animal", "human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  geom_point2(aes(subset=(node %in% unlist(cat$unfiltered_MRCAin))),color="red",size=5, shape = 25) + 
  geom_point2(aes(subset=(node %in% unlist(cat$filtered_MRCAin))),color="blue",size=4, shape = 8) +
  theme(legend.position = 'none') + geom_nodelab(aes(color =labels), size = 6)
#+ ggtitle("a. Cat") 
ggsave(cat_plot, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/cat_labs.pdf", height = 50, width = 50, limitsize = FALSE)

#groupOTU(cat$tree, focus=cat$tree$tip.label)

#dog plot subtree 4
dog = find_intro("dog", 4, c("animal", "human"))
dog_plot = ggtree(dog$tree, aes(color = labels, size = labels), layout = 'circular', branch.length = 'none') %<+% dog$all +
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  scale_size_manual(values = c(1, 0.5), labels = c("animal", "human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  geom_point2(aes(subset=(node %in% unlist(dog$unfiltered_MRCAin))),color="red",size=5, shape = 25) +
  geom_point2(aes(subset=(node %in% unlist(dog$filtered_MRCAin))),color="blue",size=4, shape = 8) +
  theme(legend.position = 'none') #+ ggtitle("b. Dog")
ggsave(dog_plot, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/dog.pdf", height = 10, width = 10)


#mink subtree 4
mink = find_intro("mink", 6, c("animal", "human"))
mink_plot = ggtree(mink$tree, aes(color = labels, size = labels), layout = 'circular', branch.length = 'none') %<+% mink$all +  
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  scale_size_manual(values = c(0.4, 0.5), labels = c("animal", "human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  geom_point2(aes(subset=(node %in% unlist(mink$unfiltered_MRCAin))),color="red",size=5, shape = 25) + 
  geom_point2(aes(subset=(node %in% unlist(mink$filtered_MRCAin))),color="blue",size=4, shape = 8)+
  #ggtitle("c. Mink") + 
  theme(legend.position = 'none') + geom_nodelab(aes(color =labels), size = 7) + geom_tiplab2(aes(color=labels))
ggsave(mink_plot, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/mink_labs.pdf", height = 350, width = 350, limitsize = FALSE)

#deer subtree 2
deer = find_intro("deer", 2, c("animal", "human"))
deer_plot = ggtree(deer$tree, aes(color = labels, size = labels), layout = 'circular', branch.length = 'none') %<+% deer$all +  
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  scale_size_manual(values = c(0.7, 0.5), labels = c("animal", "human")) +
  geom_point2(aes(subset=(node %in% unlist(deer$unfiltered_MRCAin))),color="red",size=5, shape = 25) +
  geom_point2(aes(subset=(node %in% unlist(deer$filtered_MRCAin))),color="blue",size=4, shape = 8)+
  theme(legend.position = 'none') #+ ggtitle("d. Deer")
ggsave(deer_plot, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/deer.pdf", height = 10, width = 10)

library(ggpubr)
fake_df_tree = data.frame(row1 = c("Animal-descendant Branch", "Human-descendant Branch"), row2 = 1:2)
fake_legend_tree = ggplot(fake_df_tree, aes(x=row2, y=row1)) + geom_point(aes(color = row1, shape = row1)) + 
  scale_color_manual(name = "Tree annotation", labels = c("Animal-descendant Branch", "Human-descendant Branch"), values = c("#FCA311", "#979FC4"))+
  scale_shape_manual(name = "Tree annotation",labels = c("Animal-descendant Branch", "Human-descendant Branch"), values = c('-', '-')) + 
  scale_size_manual(name = "Tree annotation",labels = c("Animal-descendant Branch", "Human-descendant Branch"), values = c(20,20)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +   theme(legend.key.size = unit(0.1, 'cm'))
as_ggplot(get_legend(fake_legend_tree))
ggsave(as_ggplot(get_legend(fake_legend_tree)), file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/leg1.pdf", device = "pdf", height = 2, width = 4)

fake_df_trans = data.frame(row1 = c("Bootstrap-filtered Transmission", "Unfiltered Transmission"), row2 = 1:2)
fake_legend_trans = ggplot(fake_df_trans, aes(x=row2, y=row1)) + geom_point(aes(color = row1, shape = row1)) + 
  scale_color_manual(name = "Transmission annotation", labels = c("Bootstrap-filtered Transmission", "Unfiltered Transmission"), values = c("blue", "red"))+
  scale_shape_manual(name = "Transmission annotation",labels = c("Bootstrap-filtered Transmission", "Unfiltered Transmission"), values = c(8, 25)) + 
  #scale_size_manual(name = "Tree annotation",labels = c("Animal-descendant Branch", "Human-descendant Branch"), values = c(20,20)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +   theme(legend.key.size = unit(0.1, 'cm'))
ggsave(as_ggplot(get_legend(fake_legend_trans)), file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/leg2.pdf", device = "pdf", height = 2, width = 4)

################ plot supp figure trees ##################
#cat plot subtree 2
#cat = find_intro("cat", 2, c("animal", "human"))
cat_plot = ggtree(cat$tree, aes(color = labels, size = labels)) %<+% cat$all +  
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  scale_size_manual(values = c(1.9, 1.5), labels = c("animal", "human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  geom_point2(aes(subset=(node %in% unlist(cat$unfiltered_MRCAin))),color="red",size=5, shape = 25) + 
  geom_point2(aes(subset=(node %in% unlist(cat$filtered_MRCAin))),color="blue",size=4, shape = 8) +
  theme(legend.position = 'none') + geom_nodelab(color = "black", size = 5) + 
  geom_tiplab(aes(color = labels), size = 7)
#+ ggtitle("a. Cat") 
ggsave(cat_plot, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/cat_2_supp.pdf", height = 450, width = 150, limitsize = FALSE)

#groupOTU(cat$tree, focus=cat$tree$tip.label)

#dog plot subtree 4
#dog = find_intro("dog", 4, c("animal", "human"))
dog_plot = ggtree(dog$tree, aes(color = labels, size = labels)) %<+% dog$all +
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  scale_size_manual(values = c(2, 1.5), labels = c("animal", "human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  geom_point2(aes(subset=(node %in% unlist(dog$unfiltered_MRCAin))),color="red",size=5, shape = 25) +
  geom_point2(aes(subset=(node %in% unlist(dog$filtered_MRCAin))),color="blue",size=4, shape = 8) +
  theme(legend.position = 'none') + geom_nodelab(color = "black", size = 5) + 
  geom_tiplab(aes(color = labels), size = 7) #+ ggtitle("b. Dog")
ggsave(dog_plot, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/dog_4_supp.pdf", height = 450, width = 150, limitsize = FALSE)


#mink subtree 6
mink = find_intro("mink", 6, c("animal", "human"))
mink_plot = ggtree(mink$tree, aes(color = labels, size = labels), branch.length = 'none', size=4) %<+% mink$all +  
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  scale_size_manual(values = c(0.4, 0.5), labels = c("animal", "human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  geom_point2(aes(subset=(node %in% unlist(mink$unfiltered_MRCAin))),color="red",size=12, shape = 25, stroke = 5) + 
  geom_point2(aes(subset=(node %in% unlist(mink$filtered_MRCAin))),color="blue",size=9, shape = 8, stroke =3)+
  #ggtitle("c. Mink") + 
  theme(legend.position = 'none') + geom_nodelab(color = "black", size = 8) + 
  geom_tiplab(aes(color = labels), size = 7)

msa = read.fasta("/Users/sana/Documents/Summer2022-mBio review/manuscript_review/mink.fasta")
chars = matrix("s", nrow = dim(msa)[1], ncol = nchar(msa$seq.text))
for (i in 1:dim(msa)[1]){
  chars[i,] = unlist(strsplit(msa$seq.text[i],""))
}
#make sure all chars are lowercase
chars = tolower(chars)
#load hits table
table = read.delim("/Users/sana/Dropbox/Summer2022-manuscript review/manuscript_review/mink_table_w_homoplasies.tsv")
table$allele_2 = tolower(table$allele_2)

mat = matrix(nrow = length(mink$tree$tip.label) ,ncol = dim(table)[1])
for(i in 1:dim(mat)[2]){
  a1=tolower(table$allele_1[i])
  a2=tolower(table$allele_2[i])
  mat[,i] = rep("Wildtype Allele", dim(mat)[1])
  ind_in_fasta = which(chars[,table$hit_position[i]] == a2)
  ids = msa$seq.name[ind_in_fasta]
  #sanity check
  if(length(ids) != length(which(ids %in% mink$tree$tip.label))){print("missinng data")
    print(i)}
  mat[match(ids, mink$tree$tip.label),i] = "Alternate Allele"
}

genotype = data.frame(mat)
rownames(genotype) = mink$tree$tip.label
colnames(genotype) = table$annotation_amino
p = gheatmap(mink_plot, genotype, width = 0.02, offset =3, colnames_position = "top", font.size = 8) +
  theme(legend.position = 'none') +
  scale_fill_manual(breaks=c("Alternate Allele", "Wildtype Allele"), 
                    values=c("#840032", "#DAD3D5"))

fake_df = data.frame(row1 = c("Animal branch/sequence", "Human branch/sequence",
                              "Bootstrap-filtered transmission", "Unfiltered transmissiom",
                              "Alternate Allele", "Wildtype Allele"), row2 = 1:6)
fake_legend = ggplot(fake_df, aes(x=row2, y=row1)) + 
  geom_point(aes(color = row1, shape = row1)) + 
  scale_color_manual(name = "Tree annotation", labels = c("Animal branch/sequence", "Human branch/sequence", "Bootstrap-filtered transmission", "Unfiltered transmissiom", "Alternate Allele", "Wildtype Allele"), values = c("orange", "#A8B2DC", "blue", "red", "darkred", "#DAD3D5")) +
  scale_shape_manual(name = "Tree annotation", labels = c("Animal branch/sequence", "Human branch/sequence", "Bootstrap-filtered transmission", "Unfiltered transmissiom", "Alternate Allele", "Wildtype Allele"), values = c(15, 15, 8, 25, 15, 15)) + 
  guides(color = guide_legend(override.aes = list(size = 25))) + 
  theme(legend.key.size = unit(2, 'cm')) + theme(legend.text=element_text(size=50)) +
  theme(legend.title=element_text(size=25))
leg1 <- get_legend(fake_legend)

p = ggarrange(p, legend.grob = leg1, legend = 'right')

ggsave(p, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/mink_6_supp.pdf", height = 1000, width = 550, limitsize = FALSE)



#deer subtree 2
deer = find_intro("deer", 2, c("animal", "human"))
deer_plot = ggtree(deer$tree, aes(color = labels, size = labels), branch.length = 'none', size = 4) %<+% deer$all +  
  scale_color_manual(values=c("#FCA311", "#A8B2DC"),labels = c("animal","human")) +
  #geom_tippoint(aes(color = trait), size = 1.5) + 
  scale_size_manual(values = c(2.7, 2.5), labels = c("animal", "human")) +
  geom_point2(aes(subset=(node %in% unlist(deer$unfiltered_MRCAin))),color="red",size=12, shape = 25, stroke = 5) +
  geom_point2(aes(subset=(node %in% unlist(deer$filtered_MRCAin))),color="blue",size=9, shape = 8, stroke = 3)+
  theme(legend.position = 'none') + geom_nodelab(color = "black", size = 10) + #theme(legend.text=element_text(size=20))
  geom_tiplab(aes(color = labels), size = 12) #+ ggtitle("d. Deer")


msa = read.fasta("/Users/sana/Documents/Summer2022-mBio review/manuscript_review/deer.fasta")
chars = matrix("s", nrow = dim(msa)[1], ncol = nchar(msa$seq.text))
for (i in 1:dim(msa)[1]){
  chars[i,] = unlist(strsplit(msa$seq.text[i],""))
}
#make sure all chars are lowercase
chars = tolower(chars)
#load hits table
table = read.delim("/Users/sana/Dropbox/Summer2022-manuscript review/manuscript_review/deer_table_w_homoplasies.tsv")
table$allele_2 = tolower(table$allele_2)

#only plot gwas hits that appear in all 10 tree replicates
mat = matrix(nrow = length(deer$tree$tip.label) ,ncol = 7)
for(i in 1:dim(mat)[2]){
  a1=tolower(table$allele_1[i])
  a2=tolower(table$allele_2[i])
  mat[,i] = rep("Wildtype Allele", dim(mat)[1])
  ind_in_fasta = which(chars[,table$hit_position[i]] == a2)
  ids = msa$seq.name[ind_in_fasta]
  #sanity check
  if(length(ids) != length(which(ids %in% deer$tree$tip.label))){print("missinng data")
    print(i)}
  mat[match(ids, deer$tree$tip.label),i] = "Alternate Allele"
}

genotype = data.frame(mat)
rownames(genotype) = deer$tree$tip.label
colnames(genotype) = table$annotation_amino[1:7]
p = gheatmap(deer_plot, genotype, width = 0.02, offset =4, colnames_position = "top", font.size = 8) +
  theme(legend.position = 'none') +
  scale_fill_manual(breaks=c("Alternate Allele", "Wildtype Allele"), 
                    values=c("#840032", "#DAD3D5")) 


p = ggarrange(p, legend.grob = leg1, legend = 'right')


ggsave(p, file = "/Users/sana/Documents/eLife_reviews/06-figure1-phylo/deer_2_supp.pdf", height = 1000, width = 550, limitsize = FALSE)

