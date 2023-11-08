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



# read input proteins and create gene_to_genome table ---------------------

faa_file <- read.fasta("vcontact_input/combined_protein_file.faa", type = "AA")
gene_to_genome_df <- data.table(protein_id = unlist(names(faa_file)),
                                contig_id = "", 
                                annot = "")

gene_to_genome_df <- gene_to_genome_df %>% 
    mutate(protein_id = str_remove(protein_id, pattern = "\\s#.*$")) %>% 
    mutate(contig_id = str_remove(protein_id, pattern= "\\_\\d*$"))


fwrite(gene_to_genome_df, file = "vcontact_input/combined_gene_to_genome_df.csv", 
       sep = ",")









