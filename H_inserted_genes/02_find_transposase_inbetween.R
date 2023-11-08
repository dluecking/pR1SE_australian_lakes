# author: dlueckin
# date: Tue Jun 27 11:59:11 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(googlesheets4)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# broder_df ---------------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(contig_category == "complete") %>% 
    filter(!is.na(manual_min_orf))


# read protein_df with annotations ----------------------------------------

prot_df <- fread("../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv")
prot_df <- prot_df %>% 
    filter(contig %in% manual_region_df$contig) %>% 
    filter(str_detect(nr_annotation, regex(".*transposase.*", ignore_case = TRUE)) |
               str_detect(interpro_annot, regex(".*transposase.*", ignore_case = TRUE)))

manual_region_df <- manual_region_df %>% 
    filter(contig %in% prot_df$contig)

manual_region_df %>% 
    select(contig, manual_min_orf, manual_max_orf)


