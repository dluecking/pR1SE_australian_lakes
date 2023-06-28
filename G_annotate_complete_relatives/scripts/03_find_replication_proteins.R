# author: dlueckin
# date: Wed Jun 21 14:37:56 2023

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


# import full df ----------------------------------------------------------

int_df <- fread("interpro_out_full.tsv")


# search replication proteins ---------------------------------------------

search_terms <- paste(c("cdc6", "orc1", "replication", "RepA", "primase",
                        "polymerase", "helicase", "jelly roll", "jelly-roll",
                        "endonuclease", " DNA-binding domain", "DNA binding domain",
                        "Origin-binding", "origin binding", "integrase"),
                      collapse = "|")
rep_df <- data.table()

for(acc in unique(int_df$contig)){
    tmp_df <- int_df %>% 
        filter(contig == acc) %>% 
        filter(grepl(pattern = search_terms, x = annotation, ignore.case = TRUE))
    rep_df <- rbind(rep_df, tmp_df)

}


fwrite(rep_df, "replication_df.tsv", sep = "\t")
