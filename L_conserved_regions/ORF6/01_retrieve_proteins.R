# author: dlueckin
# date: Wed Jul  5 15:02:02 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# get  region df ------------------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
region_df <- region_df %>% 
    filter(contig_category == "complete")


# get protein_df ----------------------------------------------------------

prot_df <- fread("../../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv") %>% 
    filter(annot == "ORF6.faa")


# loop over region_df -----------------------------------------------------

# CAREFUL
# we append, so make sure you do this only once or delete first!

for(i in 1:nrow(region_df)){
    # read protein files for this region
    seq_file <- paste0("../../G_annotate_complete_relatives/genes/", region_df$contig[i], ".fasta.faa")
    all_seqs <- read.fasta(seq_file, seqtype = "AA")
    
    # skip if no region was defined
    if(is.na(region_df$manual_min_orf[i]))
        next
    
    relevant_sequence_names <- prot_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(`orf#` >= region_df$manual_min_orf[i]) %>% 
        filter(`orf#` <= region_df$manual_max_orf[i]) %>% 
        filter(str_detect(annot, pattern = "ORF6")) %>% 
        select(id) %>% 
        pull()
    
    relevant_seqs <- all_seqs[relevant_sequence_names]
    
    # append to existing files
    out_file <- paste0("proteins/ORF6_combined.faa")
    write.fasta(sequences = relevant_seqs, 
                names = relevant_sequence_names,
                file.out = out_file,
                open = "a")
}


# only get long hits ------------------------------------------------------

all_proteins <- 
