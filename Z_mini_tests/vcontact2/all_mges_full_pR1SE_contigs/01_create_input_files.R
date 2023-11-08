# author: dlueckin
# date: Wed Jun 28 11:02:58 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# load region df ----------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives")) %>% 
    filter(contig_category == "complete") %>% 
    filter(!is.na(manual_region_length)) %>% 
    filter(integrated_state != "main_chromosome")


# copy relevant AA files --------------------------------------------------

relevant_seqs <- list()
for(i in 1:nrow(manual_region_df)){
    faa_file <- paste0("../../../FINAL_complete_relatives_proteins/aa/", manual_region_df$contig[i], ".fasta.faa")
    all_seqs <- read.fasta(faa_file)
    relevant_seqs <- append(relevant_seqs, all_seqs)
}

plasmid_seqs <- read.fasta("../plasmids/genes/archaeal_plasmids_between_20k_70k_genes.faa")
relevant_seqs <- append(relevant_seqs, plasmid_seqs)


write.fasta(relevant_seqs, names = getName(relevant_seqs),
            file = "vcontact_input/protein_file.faa")



# create gene_to_genome_file ----------------------------------------------


faa_file <- read.fasta("vcontact_input/protein_file.faa")
gene_to_genome_df <- data.table(protein_id = unlist(names(faa_file)),
                                contig_id = "", 
                                annot = "")

gene_to_genome_df <- gene_to_genome_df %>% 
    mutate(protein_id = str_remove(protein_id, pattern = "\\s#.*$")) %>% 
    mutate(contig_id = str_remove(protein_id, pattern= "\\_\\d*$"))


fwrite(gene_to_genome_df, file = "vcontact_input/gene_to_genome_df.csv", 
       sep = ",")
