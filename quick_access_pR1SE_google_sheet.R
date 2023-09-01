# author: dlueckin
# date: Wed Jul 26 14:35:54 2023

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


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))



# distance of integrase to region -----------------------------------------

region_df$integrase_distance <- NA
    
    
for(i in 1:nrow(region_df)){
    # only calc if region boundaries and an integrase are there
    if(!is.na(region_df$manual_min_orf[i]) && !is.na(region_df$integrase_orf_clean[i])){
        # check if the integrase is within the region
        if(between(region_df$integrase_orf_clean[i], 
                   left = region_df$manual_min_orf[i], 
                   right = region_df$manual_max_orf[i])){
            
            region_df$integrase_distance[i] <- 0
        }else{
            # if not inside the region, then calculate the distance
            region_df$integrase_distance[i] <- 
                min(abs(region_df$manual_max_orf[i] - region_df$integrase_orf_clean[i]), 
                    abs(region_df$manual_min_orf[i] - region_df$integrase_orf_clean[i]))
            
        }
    }
}



# add bp position of regions ----------------------------------------------


region_df$region_start_bp <- ""
region_df$region_end_bp <- ""

prot_df <- fread("G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv")


for(i in 1:nrow(region_df)){
    if(is.na(region_df$manual_max_orf[i]))
        next
    
    region_df$region_start_bp[i] <- prot_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(`orf#` == region_df$manual_min_orf[i]) %>% 
        select(start) %>% 
        pull()
    
    
    region_df$region_end_bp[i] <- prot_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(`orf#` == region_df$manual_max_orf[i]) %>% 
        select(end) %>% 
        pull()
    
}
