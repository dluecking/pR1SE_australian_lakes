# author: dlueckin
# date: Tue Oct 24 16:03:46 2023

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



# import protein cluster df -----------------------------------------------

# get only cluster 14 (which corresponds to ORF25)
cluster_df <- fread("../../F_all_pR1SE_ORFs/all_vs_all/temp_cd-hit_approach/cluster-df.tsv") %>% 
    filter(cluster == "--- Cluster 14")


all_proteins <- read.fasta("all_FINAL_COMPLETE_proteins_combined.faa")
relevant_proteins <- all_proteins[cluster_df$V1]

for(protein in relevant_proteins){
    write.fasta(protein, names = getName(protein), file.out = paste0("ORF25_proteins/", getName(protein)))
}



# do the same thing for ORF7 ----------------------------------------------

# get only cluster 14 (which corresponds to ORF25)
cluster_df <- fread("../../F_all_pR1SE_ORFs/all_vs_all/temp_cd-hit_approach/cluster-df.tsv") %>% 
    filter(cluster == "--- Cluster 16")


all_proteins <- read.fasta("all_FINAL_COMPLETE_proteins_combined.faa")
relevant_proteins <- all_proteins[cluster_df$V1]

for(protein in relevant_proteins){
    write.fasta(protein, names = getName(protein), file.out = paste0("ORF7_proteins/", getName(protein)))
}

