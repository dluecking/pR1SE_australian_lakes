# author: dlueckin
# date: Tue May 30 15:06:06 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)

# working directory -------------------------------------------------------
this_dir <- "/home/dlueckin/projects/pR1SE_australian_lakes/search_ASL/scripts"
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# read contig_df -----------------------------------------------------------

contig_df <- fread("novel_hosts_df.tsv")


# read all intput fasta files and keep relevant hits ----------------------
relevant_seqs <- list()

for(fasta_file in list.files("../tal/4_Dom/", pattern = "fasta", full.names = TRUE)){
    current_fasta <- read.fasta(fasta_file)
    relevant_seqs <- c(relevant_seqs, current_fasta[which(names(current_fasta) %in% contig_df$contig)])
}


# save to file ------------------------------------------------------------

write.fasta(relevant_seqs, names = names(relevant_seqs), file.out = "relevant_sequences.fasta")

