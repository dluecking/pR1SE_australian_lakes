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
    filter(!is.na(manual_region_length))


# read full prot_df  ------------------------------------------------------

prot_df <- fread("../../../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM.tsv")


# read input proteins and create gene_to_genome table ---------------------

gene_to_genome_df <- data.table()
for(file in list.files("../input_files/protein_files", full.names = TRUE)){
    tmp_df <- data.table(protein_id = getName(read.fasta(file)),
                         contig_id = str_remove(string = str_remove(string = file, 
                                                                    pattern = "../input_files/protein_files/"),
                                                pattern = "\\_region_only.*$"))
    tmp_df$keywords <- prot_df$annot[match(tmp_df$protein_id, prot_df$id)]
    gene_to_genome_df <- rbind(gene_to_genome_df, tmp_df)
}
fwrite(gene_to_genome_df, file = "../input_files/gene_to_genome_df.csv", 
       sep = ",")