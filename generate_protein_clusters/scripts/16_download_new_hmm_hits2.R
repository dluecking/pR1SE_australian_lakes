# author: dlueckin
# date: Wed May  3 11:18:19 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(rentrez)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# import ORF info ---------------------------------------------------------

orf_df <- fread("ORF_cutoff_table.tsv")

# import blacklist --------------------------------------------------------

blacklist <- fread("../accession_blacklist.txt", header = F) %>% 
    filter(!str_detect(V1, "#"))


# import hmmer out --------------------------------------------------------

big_df <- data.table()

for(i in 1:nrow(orf_df)){
    # read and select relevant rows
    domtbl <- fread(paste0("../08_hmmer_results/", orf_df$orf[i], "_domtblout.txt"), skip = 3, fill = TRUE) %>% 
        select(V1, V3, V4, V6, V7, V8, V16, V17)
    names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
    
    # add ORF info
    domtbl$orf <- orf_df$orf[i]
    
    # remove unuseable lines
    domtbl <- domtbl %>% 
        filter(hit != "#")
    
    # calc alignment length
    domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start
    
    # filter based on this specific orf's requirements
    domtbl$query_length <- as.numeric(domtbl$query_length)
    domtbl$score <- as.numeric(domtbl$score)
    
    domtbl <- domtbl %>% 
        filter(score >= orf_df$min_score[i]) %>% 
        filter(aln_length >= domtbl$query_length[i] * orf_df$min_alignment_length[i] / 100) %>% 
        filter(!hit %in% blacklist$V1)
    
    big_df <- rbind(big_df, domtbl)
}


# download the accessions and save to file
for(ORF in unique(big_df$orf)){
    fasta_in <- read.fasta(paste0("../12_curated_proteins/", ORF, ".faa"))
    novel_acc <- big_df %>% 
        filter(orf == ORF) %>% 
        filter(!hit %in% getName(fasta_in)) %>% 
        select(hit) %>% 
        unlist()
    
    print(paste0("current ORF: ", ORF))
    print(paste0("no of proteins before: ", length(fasta_in)))
    print(paste0("no of novel proteins: ", length(novel_acc)))


    # file out
    f_out <- paste0("../10_hmmer_proteins/", ORF, ".faa")
    # write old proteins
    #write.fasta(file.out = f_out, names = getName(fasta_in), sequences = getSequence(fasta_in))
    
    # append new records
    # if(length(novel_acc) > 0){
    #     print("downloading new accs...")
    #     novel_records <- entrez_fetch(db="protein", id=novel_acc, rettype = "fasta")
    #     write(x = novel_records, file = f_out, append = TRUE)
    # }
}



