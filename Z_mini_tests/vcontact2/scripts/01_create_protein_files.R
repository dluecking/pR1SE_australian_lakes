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


# copy relevant AA files --------------------------------------------------

relevant_seqs <- list()
for(i in 1:nrow(manual_region_df)){
    faa_file <- paste0("../../../E_check_synteny/genes/", manual_region_df$contig[i], ".fasta.faa")
    all_seqs <- read.fasta(faa_file)
    relevant_seqs <- append(relevant_seqs, all_seqs[manual_region_df$manual_min_orf[i]:manual_region_df$manual_max_orf[i]])
}

write.fasta(relevant_seqs, names = getName(relevant_seqs),
            file = "../input_files/protein_file.faa")
