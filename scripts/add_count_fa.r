#Run from base/results 
#In terminal run:
#conda activate hmsc

library(tidyverse)
library(seqinr)
library("Biostrings")

#load OTU table and fasta file
fastafile <- readDNAStringSet("~/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/base/results/DADA2_nochim.otus")
OTUtable <- read.table("~/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/base/results/DADA2_nochim.table", sep="\t", header=T, row.names=1,check.names=F)

#match seq name in OTU file with seq name in original fasta file, get sequence number and rowsums
seq_df <- OTUtable[match(row.names(OTUtable),names(fastafile)),] %>% 
  mutate(names = row.names(OTUtable), size = rowSums(OTUtable)) %>% 
  select(c(names, size))

#Add "size" and read counts to names
seq_name <- paste(seq_df$names, ";size=", seq_df$size, ";", sep = "")

#create a vector of sequences
sequence = as.list(paste(fastafile))

#write the FASTA file
write.fasta(sequences = sequence, names = seq_name, 
            as.string=FALSE, file.out="COSQ_16S.fa")

# Check seq lengths
# conda activate bioinfo
# seqmagick info base/results/COSQ_16S.fa 
# name                alignment    min_len   max_len   avg_len  num_seqs
# results/COSQ_16S.fa FALSE            160       448    254.82   1532113