# author: dlueckin
# date: Tue Sep 12 16:55:23 2023

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


WINDOW <- 0

# blast_df ----------------------------------------------------------------

blast_df <- fread("combined_blast_results.out")
names(blast_df) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", " "))


# filter ------------------------------------------------------------------

blast_df <- blast_df %>% 
    filter(pident > 95) %>% 
    filter(mismatch <= 2) %>% 
    filter(length >= 20) %>% 
    filter(qseqid != "CP019286.1") %>% 
    filter(qseqid != "AOID01000019.1")


# connect to pR1SE data sheet ---------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
region_df <- region_df %>% 
    filter(contig_category == "complete") %>% 
    filter(!is.na(region_start_bp))


# check that the spacer is within conserved region boundaries
blast_df$within_region <- FALSE
for(i in 1:nrow(blast_df)){
    spacer_position <- blast_df$qstart[i]
    contig <-  blast_df$qseqid[i]
    region_start <-  region_df$region_start_bp[region_df$contig == contig]
    region_end <- region_df$region_end_bp[region_df$contig == contig]
    if(between(x = spacer_position, left = region_start - WINDOW, right = region_end + WINDOW)){
        blast_df$within_region[i] <- TRUE
    }
}

true_hits <- blast_df %>% filter(within_region)

region_df$CIRSPR_targeted_by <- ""
for(i in 1:nrow(region_df)){
    hosts <- blast_df %>% 
        filter(qseqid == region_df$contig[i]) %>% 
        filter(within_region) %>% 
        select(host_species) %>% 
        unique() %>% 
        unlist()
    hosts <- paste(hosts, collapse = "; ")
    
    region_df$CIRSPR_targeted_by[i] <- hosts
    
    
}
write(region_df$CIRSPR_targeted_by, "CRISPR_targeted_by5k.tsv")
