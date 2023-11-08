# author: dlueckin
# date: Wed Oct 11 11:36:09 2023

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

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))


# load blast df -----------------------------------------------------------

blast_df <- fread("vcontact2_out/merged.self-diamond.tab")
names(blast_df) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", " "))

blast_df$q_contig_id <- str_remove(blast_df$qseqid, "\\_\\d.*$")
blast_df$s_contig_id <- str_remove(blast_df$sseqid, "\\_\\d.*$")


# explore -----------------------------------------------------------------
# it is complicated: I need to find the proteins shared between a couple of connections:

# CP035122.1 - AOIO01000049.1
tmp_df <- blast_df %>% 
    filter(q_contig_id == "CP035122.1") %>% 
    filter(s_contig_id == "AOIO01000049.1")

# CP128375.1 - CP031299.1
# CP128375.1 - JDTG01000019.1
connect_to_5 <- blast_df %>% 
    filter(q_contig_id == "CP128375.1") %>% 
    filter(s_contig_id == "CP031299.1" | s_contig_id == "JDTG01000019.1")


# CP129993.1 - CP019286.1

# CP129993.1 - CP019329.1

# infact, there are more connections...but I just want to see if the connections 
# hit proteins in the core region of the second contig









