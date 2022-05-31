#########
#cleanup animal files download from gisaid


library(phylotools)

msa1 = read.fasta("1646449641959.sequences.fasta")
metadata1 = read.delim("1646449641959.metadata.tsv")
dim(msa1)
dim(metadata1)



msa2 = read.fasta("1646449907104.sequences.fasta")
metadata2 = read.delim("1646449907104.metadata.tsv")

dim(msa2)
dim(metadata2)


#bind them:
all_msa = rbind(msa1,msa2)
all_metadata = rbind(metadata1,metadata2)

# remove dupes:
dupes = which(duplicated(all_msa$seq.name))
if(length(dupes)>0){all_msa = all_msa[-dupes,]}

dupes = which(duplicated(all_metadata[,1]))
if(length(dupes)>0){all_metadata = all_metadata[-dupes,]}


#keep only the metadata of existing sequences:
all_metadata = all_metadata[which(all_metadata[,1] %in% all_msa$seq.name),]

msa = all_msa
metadata = all_metadata


# remove unwanted hosts
hosts = unique(metadata$host)
to_remove = c("unknown","Environment")
metadata = metadata[-which(metadata$host %in% to_remove),]
msa = msa[which(msa$seq.name %in% metadata[,1]),]


library(stringr)
# remove incomplete dates
metadata = metadata[which(str_count(metadata$date,"-") == 2),]
msa = msa[which(msa$seq.name %in% metadata[,1]),]



# remove too many n's
msa = msa[-which(str_count(msa$seq.text,"N")>500),]
metadata = metadata[which(metadata[,1] %in% msa$seq.name),]


#remove dashes
msa$seq.text = str_remove_all(msa$seq.text,"-")

library(seqinr)
#write fasta
write.fasta(names = as.list(msa$seq.name), sequences = as.list(msa$seq.text), file.out = "animal_msa_gisaid_04032022.fasta")
write.table(metadata, file = "animal_metadata_gisaid_04032022.tsv", row.names = FALSE, sep = "\t")


