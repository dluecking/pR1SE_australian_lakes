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



# add direction of region -------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
region_df <- region_df %>% 
    filter(contig_category == "complete")


prot_df <- fread("G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv") %>% 
    filter(annot == "ORF6.faa" | annot == "ORF24.faa")


for(i in 1:nrow(region_df)){
    orf24_pos <- prot_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(annot == "ORF24.faa")

    orf6_pos <- prot_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(annot == "ORF6.faa")
    
    if(nrow(orf24_pos) > 1 | nrow(orf6_pos) > 1){
        
        print("weird")
        print(i)
        
        orf24_pos <- prot_df %>% 
            filter(contig == region_df$contig[i]) %>% 
            filter(annot == "ORF24.faa") %>% 
            sample_n(1)
        
        orf6_pos <- prot_df %>% 
            filter(contig == region_df$contig[i]) %>% 
            filter(annot == "ORF6.faa") %>% 
            sample_n(1)
    }
    region_df$region_direction[i] <- ifelse(orf24_pos$start[1] >= orf6_pos$start[1],
                                         "forward",
                                         "reverse")
    
}
# write(region_df$region_direction, "quui.txt")
