library(stringr)
#read mink hit table
mink_hits = read.delim("mink_table_w_homoplasies.tsv")
deer_hits = read.delim("deer_table_w_homoplasies.tsv")

### all positions are in alignment numbering, NOT REFERENCE 
what_base = function(species, gwas_hit, wildtype, alternate){
  df = data.frame(species = vector(),
                  hit_pos_aln = character(),
                  replicate = numeric(),
                  out_node_number = numeric(),
                  in_node_number = numeric(),
                  site_mutated_in_intro = character(),
                  parent_allele = character(),
                  child_allele = character(),
                  same_sub_as_gwas = character())
  for(i in 1:10){
    #read intro sequences table
    tab = read.delim(paste0(species, "_", i, "_intro_sequences_table_human2animal.tsv"))
    print("replicate:")
    print(i)
    
    
    for(j in 1:dim(tab)[1]){
      
      #if there are no differences skip to next intro
      if(as.numeric(tab$ndiff[j]) == 0){
        df[nrow(df) + 1, ] = c(species, 
                               gwas_hit,
                               i,
                               tab$out_node_number[j],
                               tab$in_node_number[j],
                               "no",
                               tolower(substr(tab$out_sequence[j], gwas_hit, gwas_hit)),
                               tolower(substr(tab$in_sequence[j], gwas_hit, gwas_hit)),
                               "-")
        next}
      sites = unlist(str_split(tab$sites[j], "/"))
      
      #check if the hit is in the list of segregating sites at all
      if(!as.character(gwas_hit) %in% sites){
        df[nrow(df) + 1, ] = c(species, 
                               gwas_hit,
                               i,
                               tab$out_node_number[j],
                               tab$in_node_number[j],
                               "no",
                               tolower(substr(tab$out_sequence[j], gwas_hit, gwas_hit)),
                               tolower(substr(tab$in_sequence[j], gwas_hit, gwas_hit)),
                               "-")
        next
      }
      if(as.character(gwas_hit) %in% sites){
        flag = "no"
        if(tolower(substr(tab$out_sequence[j], gwas_hit, gwas_hit)) == wildtype){
          if(tolower(substr(tab$in_sequence[j], gwas_hit, gwas_hit)) == alternate){
            flag = "yes"
          }
        }
        df[nrow(df) + 1, ] = c(species, 
                               gwas_hit,
                               i,
                               tab$out_node_number[j],
                               tab$in_node_number[j],
                               hit_in_intro = "yes",
                               tolower(substr(tab$out_sequence[j], gwas_hit, gwas_hit)),
                               tolower(substr(tab$in_sequence[j], gwas_hit, gwas_hit)),
                               flag)
      }
    }
  }
  #write output
  write.table(df, file = paste0(species, "_hit_", gwas_hit, "_summary.tsv"), row.names = FALSE, sep = "\t")
}

## sanity checks:
typeof(mink_hits$allele_1)
typeof(mink_hits$allele_2)
typeof(deer_hits$allele_1)
typeof(mink_hits$allele_2)

mapply(what_base, "mink", as.numeric(mink_hits$hit_position), tolower(mink_hits$allele_1), tolower(mink_hits$allele_2))
mapply(what_base, "deer", as.numeric(deer_hits$hit_position), tolower(deer_hits$allele_1), tolower(deer_hits$allele_2))

