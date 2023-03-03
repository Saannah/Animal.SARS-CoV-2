### plot intro mutations distribution
setwd("/users/sana/Documents/eLife_reviews/08-mutation_spectrum_intros/")
## load gwas significant hits
mink_hits = read.delim("/Users/sana/Documents/AIM1&2/manuscript/tables/mink_table_w_homoplasies.tsv")
deer_hits = read.delim("/Users/sana/Documents/AIM1&2/manuscript/tables/deer_table_w_homoplasies.tsv")

fill_gaps = function(df, species){
  #convert positions to ref
  df$pos_in_ref = mapping[df$site]
  #remove indels
  if(length(which(df$pos_in_ref == -1))>0){
    df = df[-which(df$pos_in_ref == -1),]
  }
  df$species = species
  # #length of raw ref is 29891
  return(df)
}

################## animal to human direction ####################
cat = fill_gaps(read.delim("cat_mutation_summary_animal2human.tsv"), "cat")
dog = fill_gaps(read.delim("dog_mutation_summary_animal2human.tsv"), "dog")
mink = fill_gaps(read.delim("mink_mutation_summary_animal2human.tsv"), "mink")
deer = fill_gaps(read.delim("deer_mutation_summary_animal2human.tsv"), "deer")

all = rbind(fill_gaps(cat, "cat"), fill_gaps(dog, "dog"),fill_gaps(mink, "mink"), fill_gaps(deer, "deer"))

# cat$pos_in_ref = mapping[cat$site]
# if(length(which(cat$pos_in_ref == -1))>0){
#   cat = cat[-which(cat$pos_in_ref == -1),]
# }



library(ggplot2)
p1 = ggplot(cat, aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 300)) + 
  ggtitle(paste0("Cat-Site Frequency Spectrum (Animal to Human direction)")) + ylab("Frequency") + xlab("Mutation Count")
ggsave(p1, file = "cat-sfs_animal2human.pdf", device = "pdf", height = 5, width = 15)

p2 = ggplot(dog, aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 300)) + 
  ggtitle(paste0("Dog-Site Frequency Spectrum (Animal to Human direction)")) + ylab("Frequency") + xlab("Mutation Count")
ggsave(p2, file = "dog-sfs_animal2human.pdf", device = "pdf", height = 5, width = 15)

p3 = ggplot(all[which(all$species == "mink"), ], aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 300)) + 
  ggtitle(paste0("Mink-Site Frequency Spectrum (Animal to Human direction)")) + ylab("Frequency") + xlab("Mutation Count")
ggsave(p3, file = "mink-sfs_animal2human.pdf", device = "pdf", height = 5, width = 15)

p4 = ggplot(all[which(all$species == "deer"), ], aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 300)) + 
  ggtitle(paste0("Deer-Site Frequency Spectrum (Animal to Human direction)")) + ylab("Frequency") + xlab("Mutation Count")
ggsave(p4, file = "deer-sfs_animal2human.pdf", device = "pdf", height = 5, width = 15)

#ggplot(all, aes(x=pos_in_ref, y=species, fill = count)) + geom_tile()

################## human to animal direction ####################
#load human to animal data
cat = fill_gaps(read.delim("cat_mutation_summary_human2animal.tsv"), "cat")
dog = fill_gaps(read.delim("dog_mutation_summary_human2animal.tsv"), "dog")
mink = fill_gaps(read.delim("mink_mutation_summary_human2animal.tsv"), "mink")
deer = fill_gaps(read.delim("deer_mutation_summary_human2animal.tsv"), "deer")

all = rbind(cat, dog, mink, deer)


library(ggplot2)
p1 = ggplot(cat, aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 600)) + 
  ggtitle(paste0("Cat-Site Frequency Spectrum (Human to Animal direction)")) + ylab("Frequency") + xlab("Mutation Count")
ggsave(p1, file = "cat-sfs_human2animal.pdf", device = "pdf", height = 5, width = 15)

p2 = ggplot(dog, aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 300)) + 
  ggtitle(paste0("Dog-Site Frequency Spectrum (Human to Animal direction)")) + ylab("Frequency") + xlab("Mutation Count")
ggsave(p2, file = "dog-sfs_human2animal.pdf", device = "pdf", height = 5, width = 15)


p3 = ggplot(mink, aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 500)) + 
  ggtitle(paste0("Mink-Site Frequency Spectrum (Human to Animal direction)")) + 
  ylab("Frequency") + xlab("Mutation Count") + theme_bw() +
  geom_vline(xintercept=47, color = 'orange', linetype=2) + geom_text(x=47, y=0.05, label = "L219V",color = "darkblue", angle=90, size=2)+
  geom_vline(xintercept=15, color = 'orange', linetype=2) + geom_text(x=15, y=0.05, label = "N501T",color = "darkblue", angle = 90, size=2)+
  geom_vline(xintercept = 51, color = 'orange', linetype = 2) + geom_text(x=51, y=0.05, label = "G37E",color = "darkblue", angle =90,size=2)
  
  

ggsave(p3, file = "mink-sfs_human2animal.pdf", device = "pdf", height = 5, width = 15)


p4 = ggplot(deer, aes(x=count)) + geom_density(alpha = 0.4, fill = "blue") + xlim(c(-1, 300)) + 
  ggtitle(paste0("Deer-Site Frequency Spectrum (Human to Animal direction)")) + ylab("Frequency") + xlab("Mutation Count") + 
  geom_vline(data=deer[which(deer$pos_in_ref %in% deer_hits$pos_in_ref),], aes(xintercept = count), color = 'orange', linetype ='twodash')
ggsave(p4, file = "deer-sfs_human2animal.pdf", device = "pdf", height = 3, width = 15)



vec = vector()
for(i in 1:dim(deer_msa)[1]){
  vec = c(vec, substr(deer_msa$seq.text[i], 210, 210))
}
deer_msa$seq.name[which(vec == "t")]

# ################ test heatmap ##############
library(plotly)
p=ggplot(all, aes(x=pos_in_ref, y=species, fill = count)) + geom_tile() +
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  theme_bw() +
 scale_fill_gradient(low="white",high="darkblue")
ggsave(p, file ="heatmap.pdf", device = "pdf", height = 20, width = 3000, limitsize = FALSE)
