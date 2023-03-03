##script to compare the parent end of a transmission with its child and report the mutations

library(stringr)

for(species in c("cat", "dog", "mink", "deer")){
  
  variant_sites = data.frame(site = character(),
                             count = numeric())
  for(i in 1:10){
    tab = read.delim(paste0(species, "_", i, "_intro_sequences_table_human2animal.tsv"))
    print("replicate:")
    print(i)
    #parse variant sites
    for(j in 1:dim(tab)[1]){
      #print("intro:")
      #print(j)
      if(tab$ndiff[j] == 0){next}
      sites = unlist(str_split(tab$sites[j], "/"))
      present_sites = sites[which(sites %in% variant_sites$site)]
      absent_sites = sites[which(!sites %in% variant_sites$site)]
      
      if(length(present_sites) > 0){
        variant_sites$count[which(variant_sites$site %in% present_sites)] = variant_sites$count[which(variant_sites$site %in% present_sites)] + 1
        #variant_sites$reps[which(variant_sites$site %in% present_sites)] = paste0(variant_sites$reps[which(variant_sites$site %in% present_sites)], i, collapse = "/")
      }
      if(length(absent_sites) > 0){
        variant_sites = rbind(variant_sites,
                              data.frame(site = absent_sites, count = 1))
      }
      if(length(which(variant_sites$site == ""))>0){print(j)}
    }
  }
  
  write.table(variant_sites, file = paste0(species, "_mutation_summary_human2animal.tsv"), row.names = FALSE, sep = "\t", quote =TRUE)
}
