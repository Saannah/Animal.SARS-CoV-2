setwd("/Users/sana/Documents/eLife_reviews/09-permutation_on_new_trees/animal2human/")

plot_permutations = function(direction){
  summary = data.frame(species = vector(),
                  observed_sd = vector(),
                  pooled_permutated_sd = vector(),
                  average_sd = vector())
  
  for(species in c("deer", "cat", "dog", "mink")){
    mink_counts = read.delim(paste0("/Users/sana/Documents/eLife_reviews/00-intros_corrected/counts_", direction, ".tsv"))
    mink_counts = mink_counts[which(mink_counts$species == species),]

    mink_permutations = data.frame(species = species, subtree_number = 1:10, intro_count_unfiltered = mink_counts$unfiltered_count)
    
    num = vector()
    color = vector()
    draws = vector()
    sd = vector()
    for (i in 1:10){
      p = read.delim(paste0(species, "_", i, "permutations.tsv"), header = FALSE)
      #pvalue = length(which(p[,2] > mink_permutations$intro_count_filtered[i]))/1000
      #print(pvalue)
      p_sd = sd(p[,1])
      sd = c(sd, p_sd)
      draws = c(draws ,c(p[,1], mink_permutations$intro_count_unfiltered[i]))
      num = c(num, rep(i, 1001))
      color = c(color, c(rep(1, 1000), 2))
    }
    print("observed sd:")
    print(sd(mink_permutations$intro_count_unfiltered))
    print("average sd:")
    print(mean(sd))
    
    df = data.frame(draws = draws, num = num, color = color)
    
    print("pooled sd")
    print(sd(df$draws[which(df$color == 1)]))
    
    plot = ggplot(df, aes(x = num, y = draws, fill = color)) + geom_point(aes(color = color)) + ggtitle(species)
    ggsave(plot, file = paste0(species, "_all_transision_permutations.pdf"), device = "pdf")
    summary[nrow(summary) + 1, ] = c(species, 
                                     sd(mink_permutations$intro_count_unfiltered),
                                     sd(df$draws[which(df$color == 1)]),
                                     mean(sd))
  }
  write.table(summary, file = "permutation_summary.tsv", sep = "\t", row.names = FALSE)
  }

plot_permutations("animal_2_human")


setwd("/users/sana/Documents/eLife_reviews/09-permutation_on_new_trees/human2animal/")
plot_permutations("human_2_animal")



### some checks 
## check observed sd:

for(species in c("cat", "dog", "mink", "deer")){
  mink_counts = read.delim("/Users/sana/Documents/eLife_reviews/00-intros_corrected/counts_animal_2_human.tsv")
  mink_counts = mink_counts[which(mink_counts$species == species),]
  mink_permutations = data.frame(species = species, subtree_number = 1:10, intro_count_unfiltered = mink_counts$unfiltered_count)
  print(species)
  print("observed sd")
  print(sd(mink_permutations$intro_count_unfiltered))
}
