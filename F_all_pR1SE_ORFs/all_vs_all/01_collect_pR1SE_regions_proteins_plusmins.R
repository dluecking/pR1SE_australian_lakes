# author: dlueckin
# date: Tue Aug  8 15:25:06 2023

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


# load df -----------------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(contig_category == "complete") %>% 
    filter(!is.na(manual_min_orf))


# collect region proteins +/- 10 ------------------------------------------

for(i in 1:nrow(manual_region_df)){
    fasta_file <- paste0("../../FINAL_complete_relatives_proteins/aa/", manual_region_df$contig[i], ".fasta.faa")
    all_proteins <- read.fasta(fasta_file)
    
    start <- max(0, manual_region_df$manual_min_orf[i]-10)
    end <- min(manual_region_df$manual_max_orf[i]+10, length(all_proteins))
    
    relevant_proteins <- all_proteins[start:end]
    write.fasta(relevant_proteins, getName(relevant_proteins), file = paste0("proteins/",  manual_region_df$contig[i], "_relevant.faa"))
    
}
