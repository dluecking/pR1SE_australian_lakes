# author: dlueckin
# date: Mon Aug 14 13:07:06 2023

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


# import prot_df ----------------------------------------------------------

prot_df <- fread("prot_df_with_acc_and_id.csv")


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))


# filter, removing the known hosts ----------------------------------------

prot_df <- prot_df %>% 
    filter(!str_detect(string = organism_acc,
                       pattern = paste(str_remove(
                           string = region_df$contig, 
                           pattern = "\\.\\d*$"), 
                           collapse = "|"))) %>% 
    distinct() 


# remove one exception, 2 we already found but discarded
prot_df <- prot_df %>% 
    filter(organism_acc != "NZ_CP050018",
           organism_acc != "NZ_JAANTI010000001",
           organism_acc != "NZ_PSYY01000019")

# remove anythin with NA as organism, those are ones we most likley detected already
prot_df <- prot_df %>% 
    filter(organism_acc != "")

fwrite(prot_df, "filtered_prot_df.tsv")