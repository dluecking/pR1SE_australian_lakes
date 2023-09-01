# author: dlueckin
# date: Wed Jun 21 14:37:56 2023

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


# import full df ----------------------------------------------------------

int_df <- fread("interpro_out_full.tsv")


# search replication proteins ---------------------------------------------

search_terms <- paste(c("cdc6", "orc1", "replication", "RepA", "primase",
                        "polymerase", "helicase", "jelly roll", "jelly-roll",
                        "endonuclease", " DNA-binding domain", "DNA binding domain",
                        "Origin-binding", "origin binding", "PCNA",
                        "sliding clamp", "clamp", "RadA", "recombinase", "WhiP",
                        "GTPase", "rfc"),
                      collapse = "|")
rep_df <- data.table()

for(acc in unique(int_df$contig)){
    tmp_df <- int_df %>% 
        filter(contig == acc) %>% 
        filter(grepl(pattern = search_terms, x = annotation, ignore.case = TRUE))
    rep_df <- rbind(rep_df, tmp_df)

}

fwrite(rep_df, "replication_df.tsv", sep = "\t")


# subset df by region -----------------------------------------------------

# plus minus
SEARCH_REGION <- 20


manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))

rep_df <- left_join(rep_df, manual_region_df %>% select(contig, manual_min_orf, manual_max_orf))
# remove nas
rep_df <- rep_df %>% 
    filter(!is.na(manual_max_orf))
# keep only region +/- 20
subset_df <- rep_df[rep_df$orf >= rep_df$manual_min_orf - SEARCH_REGION & rep_df$orf <= rep_df$manual_max_orf + SEARCH_REGION]
# remove the winged helix
subset_df <- subset_df %>% filter(!str_detect(annotation, pattern = "Winged helix"))

fwrite(subset_df, "replication_subset_df.tsv", sep = "\t")
