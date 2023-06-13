# author: dlueckin
# date: Wed Jun  7 17:23:50 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(dplyr)
library(rentrez)

# local
# setwd("/run/user/1000/gvfs/sftp:host=linux-desktop-1.mpi-bremen.de,user=dlueckin/home/dlueckin/projects/pR1SE_australian_lakes/A_generate_protein_clusters/scripts")

# CONSTANTS AND PATHS -----------------------------------------------------

ORF_TABLE_PATH <- "../../helper_files/ORF_cutoff_table.csv"
BLACKLIST_PATH <- "../accession_blacklist.txt"


# user input --------------------------------------------------------------

options <- commandArgs(trailingOnly = TRUE)
ORF <- options[1]
DOMTBLE_FILE <- options[2]
ORIGINAL_FASTA_FILE <- options[3]
OUTPUT_FILE <- options[4]

# ORF <- "ORF10"
# DOMTBLE_FILE <- "tmp_ORF10/ORF10_hmmsearch_result.domtableout.txt"
# ORIGINAL_FASTA_FILE <- "../03_blastp_proteins/ORF10_blastp_proteins.faa"
# OUTPUT_FILE <- "test.faa"

# import ORF info ---------------------------------------------------------

orf_df <- fread(ORF_TABLE_PATH)
ORF_SCORE_CUTOFF <- orf_df$min_score[orf_df$orf == ORF]
ORF_ALN_LEN_CUTOFF <- orf_df$min_alignment_length[orf_df$orf == ORF]
rm(orf_df)


# import blacklist --------------------------------------------------------

blacklist <- fread(BLACKLIST_PATH, header = F) %>% 
    filter(!str_detect(V1, "#"))


# read domtble file -------------------------------------------------------
# read and select relevant rows
domtbl <- fread(DOMTBLE_FILE, skip = 3, fill = TRUE) %>% 
    select(V1, V3, V4, V6, V7, V8, V16, V17)
names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")

# remove unuseable lines
domtbl <- domtbl %>% 
    filter(hit != "#")

# calc alignment length
domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start

# filter based on this specific orf's requirements
domtbl$query_length <- as.numeric(domtbl$query_length)
domtbl$score <- as.numeric(domtbl$score)

domtbl <- domtbl %>% 
    filter(score >= ORF_SCORE_CUTOFF) %>% 
    filter(aln_length >= domtbl$query_length[1] * ORF_ALN_LEN_CUTOFF / 100) %>% 
    filter(!hit %in% blacklist$V1)

# remove based on the old protein file
fasta_in <- read.fasta(ORIGINAL_FASTA_FILE)
novel_acc <- domtbl %>%
    filter(!hit %in% getName(fasta_in)) %>%
    select(hit) %>%
    unlist()

print(paste0("current ORF: ", ORF))
print(paste0("no of proteins before: ", length(fasta_in)))
print(paste0("no of novel proteins: ", length(novel_acc)))
print(novel_acc)

# write old proteins
write.fasta(file.out = OUTPUT_FILE, names = getName(fasta_in), sequences = getSequence(fasta_in))

# append new records
if(length(novel_acc) > 0){
    print("downloading new accs...")
    novel_records <- entrez_fetch(db="protein", id=novel_acc, rettype = "fasta")
    write(x = novel_records, file = OUTPUT_FILE, append = TRUE)
}
