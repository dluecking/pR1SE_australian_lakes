# author: dlueckin
# date: Wed May 31 15:23:12 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gggenes)
library(googlesheets4)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# create automated region borders df --------------------------------------

gene_df <- fread("gene_df.tsv")

auto_df <- gene_df %>% 
    group_by(contig) %>% 
    summarise(min_orf = min(`orf#`), max_orf = max(`orf#`))
auto_df$region_length <- auto_df$max_orf - auto_df$min_orf
auto_df$manual_min_orf <- NA
auto_df$manual_max_orf <- NA
auto_df$manual_region_length <- NA

fwrite(auto_df, "automated_region_borders.tsv")

# I ASSUME YOU HAVE SETUP THE MANUAL_REGION_DF BASED ON THE AUTOMATED ONE BUT WITH MANUAL CURATION
rm(auto_df)



# import df and filter  --------------------------------------------------


manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))


# categorize contigs ------------------------------------------------------

manual_region_df$number_of_unique_hits  <- 0
for(i in 1:nrow(manual_region_df)){
    manual_region_df$number_of_unique_hits[i] <-  gene_df %>% 
        filter(contig == manual_region_df$contig[i]) %>% 
        summarise(length(unique(annot))) %>% 
        pull() - 1
}

manual_region_df$contig_category <- "fragment"
manual_region_df$contig_category[manual_region_df$number_of_unique_hits >= 6] <- "complete"



# save complete ones to folder --------------------------------------------

complete_relatives <- manual_region_df %>% 
    filter(contig_category == "complete")

for(acc in unique(complete_relatives$contig)){
    cmd = paste0("cp ../contigs/", acc, ".fasta", 
                 " ../../FINAL_complete_relatives/",  acc, ".fasta")
    system(cmd)
}
