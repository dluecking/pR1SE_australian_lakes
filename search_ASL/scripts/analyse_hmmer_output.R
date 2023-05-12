# author: dlueckin
# date: Tue May  9 13:14:32 2023

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



# explore -----------------------------------------------------------------

lakes <- fread("list_of_lakes.txt", header = FALSE)
orf_df <- fread("ORF_cutoff_table.tsv")



big_df <- data.table()

for(file in list.files("../01_hmm_results", full.names = TRUE, pattern = "domtblout.txt")){
    is_not_empty <- as.numeric(str_split(system(command = paste0("wc -l ", file), intern = TRUE), " ")[[1]][1]) > 13
    
    if(is_not_empty){
        domtbl <- fread(file, skip = 3, fill = TRUE) %>% 
            select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
            filter(V1 != "#")
        
        names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
        
        current_lake <- str_extract(file, pattern = paste(lakes$V1, collapse = "|"))
        domtbl$lake <- current_lake
        current_ORF <- str_extract(file, pattern = "ORF[^\\_]*")
        
        # calc alignment length
        domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start
        
        # prepare data for later filtering
        domtbl$query_length <- as.numeric(domtbl$query_length)
        domtbl$score <- as.numeric(domtbl$score)
        
        domtbl <- domtbl %>% 
            filter(score >= 50)
            # filter(score >= orf_df[orf_df$orf == current_ORF]$min_score) %>% 
            # filter(aln_length >= orf_df[orf_df$orf == current_ORF]$min_alignment_length)
        
        big_df <- rbind(big_df, domtbl)
    }
}

big_df$contig <- str_remove(big_df$hit, pattern = "\\_length.*")
contig_df <- data.table(table(big_df$contig))
