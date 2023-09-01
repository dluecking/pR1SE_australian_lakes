# author: dlueckin
# date: Mon Aug 14 15:09:51 2023

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

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))


manual_region_df$distance <- manual_region_df$max_orf - manual_region_df$min_orf







# I need the contigs, for which the hits are spread out -------------------

gene_df <- fread("gene_df.tsv")


contig_df <- data.table(contig = unique(gene_df$contig),
                        min_hit = 0)
contig_min_hit <- 0
contig_df$max_hit <- 1
for(i in 1:nrow(contig_df)){
    
    
    
    
}