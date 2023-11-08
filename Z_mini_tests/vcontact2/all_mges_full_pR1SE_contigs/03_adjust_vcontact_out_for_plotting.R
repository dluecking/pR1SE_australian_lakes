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

genome <- fread("vcontact2_out/genome_by_genome_overview.csv")
genome$label <- "Archaeal virus"

# all plasmids get label = plasmids
plasmid_ids <- fread("../plasmids/plasmid_accessions_to_download.txt", header = F)
genome$label[genome$Genome %in% plasmid_ids$V1] <- "Archaeal plasmid"


# all pR1SE elements get this title
manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
genome$label[genome$Genome %in% manual_region_df$contig] <- "pR1SE plasmid"

table(genome$label)


genome$length <- manual_region_df$sequence_length[match(genome$Genome, manual_region_df$contig)]
genome$length[is.na(genome$length)] <- 0
fwrite(genome, "genome_by_genome_overview.csv", col.names = F, sep = ",")
