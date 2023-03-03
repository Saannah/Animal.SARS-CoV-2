
######## parse gwas in intros tables
mink_df = data.frame(species = vector(),
                     hit_pos = vector(),
                     pos_in_ref = vector(),
                     site_appearance_count = vector(),
                     gwas_sub_appearance_count = vector())

for(gwas_hit in mink_hits$hit_position){
  #load table
  tab = read.delim(paste0("mink_hit_", gwas_hit, "_summary.tsv"))
  print("number of times site is mutated:")
  print(length(which(tab$site_mutated_in_intro == "yes")))
  print("number of times its exactly the gwas hit:")
  print(length(which(tab$same_sub_as_gwas == "yes")))
  mink_df[nrow(mink_df)+1, ] = c("mink",
                                 gwas_hit,
                                 mapping[gwas_hit],
                                 length(which(tab$site_mutated_in_intro == "yes")),
                                 length(which(tab$same_sub_as_gwas == "yes"))) 
  #sanity check
  if(length(which(tab$site_mutated_in_intro=="yes")) == mink$count[which(mink$site == gwas_hit)]){
    print("sanity check is ok")
  } else{
    print("ERROR")
  }
  
}


deer_df = data.frame(species = vector(),
                     hit_pos = vector(),
                     pos_in_ref = vector(),
                     site_appearance_count = vector(),
                     gwas_sub_appearance_count = vector())

for(gwas_hit in deer_hits$hit_position){
  #load table
  tab = read.delim(paste0("deer_hit_", gwas_hit, "_summary.tsv"))
  print("number of times site is mutated:")
  print(length(which(tab$site_mutated_in_intro == "yes")))
  
  deer_df[nrow(deer_df)+1, ] = c("deer",
                                 gwas_hit,
                                 mapping[gwas_hit],
                                 length(which(tab$site_mutated_in_intro == "yes")),
                                 length(which(tab$same_sub_as_gwas == "yes"))) 
  
  if(length(which(tab$site_mutated_in_intro == "yes"))==0){
    print("hit not in intros")
    next
  }
  print("number of times its exactly the gwas hit:")
  print(length(which(tab$same_sub_as_gwas == "yes")))
  #sanity check
  if(length(which(tab$site_mutated_in_intro=="yes")) == deer$count[which(deer$site == gwas_hit)]){
    print("sanity check is ok")
  } else{
    print("ERROR")
  }
}

write.table(mink_df, file = "mink_gwas_intro_report.tsv", sep ="\t", row.names = FALSE)
write.table(deer_df, file = "deer_gwas_intro_report.tsv", sep ="\t", row.names = FALSE)


