# author: dlueckin
# date: Fri Sep 22 14:43:44 2023

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


# open ntw file -----------------------------------------------------------

ntw <- fread("vcontact2_out/c1.ntw")
fwrite(ntw, "c1.csv", sep = ",", col.names = F)



# open genome file --------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))

genome <- fread("vcontact2_out/genome_by_genome_overview.csv")
genome <- genome %>% 
    mutate("order" = manual_region_df$order[match(Genome, manual_region_df$contig)]) %>% 
    mutate("replication_group" = manual_region_df$replication_group[match(Genome, manual_region_df$contig)]) 
fwrite(genome, "genome_by_genome_overview.csv", col.names = F, sep = ",")
