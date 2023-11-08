# author: dlueckin
# date: Fri Sep 22 17:12:01 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(googlesheets4)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))

# region_df
region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))



# explore -----------------------------------------------------------------

# NCBI search term:
# "archaeal"[Organism] OR "Archaea"[Organism]) AND "plasmid"[Filter] NOT "shotgun"[All Fields] 

# Read the file into a character vector
file_lines <- readLines("nuccore_result.txt")

# Create the data frame
plasmid_df <- data.table(
    name = as.character(),
    length = as.numeric(),
    circularity = as.character(),
    accession = as.character()
)


for(block in 1:(length(file_lines)/4)){
    line_blank <- block * 4 - 3
    name <- str_remove(file_lines[block * 4 - 2], "\\d*\\.\\s")
    length <- as.numeric(str_remove(string = str_extract(string = file_lines[block * 4 - 1],
                                                         pattern = "[\\d,]*"),
                         pattern = ","))

    circularity <- str_detect(file_lines[block * 4 - 1], pattern = "circular")
    
    accession <- str_remove(string = file_lines[block * 4 - 0], pattern = "\\s.*$")
    
    tmp_df <- data.table(
        name = name,
        length = length,
        circularity = circularity,
        accession = accession
    )
    plasmid_df <- rbind(plasmid_df, tmp_df)
}


# remove viruses, non plasmid containing, get hostname, remove duplicates based on length and host
# remove duplicates after NZ_ removal
# check if its in pR1SE df
# remove pR1SE contigs
cleaned_df  <- plasmid_df %>% 
    filter(str_detect(plasmid_df$name, "plasmid")) %>% 
    mutate(host = str_remove(name, "\\splasmid.*$")) %>% 
    filter(!str_detect(name, "virus")) %>% 
    filter(!str_detect(name, "Virus")) %>% 
    mutate(accession = str_remove(accession, "^NZ\\_")) %>% 
    distinct() %>% 
    mutate(is_pR1SE = accession %in% region_df$contig) %>% 
    filter(!is_pR1SE) %>% 
    filter(!str_detect(name, "shotgun"))

a <- cleaned_df %>% 
    mutate(host_family = str_extract(host, "^[A-z]*"))
table(a$host_family)


ggplot(a, aes(x = "plasmid sizes", y = length)) +
    geom_violin() +
    geom_jitter() +
    theme_cowplot() 
    


accessions <- cleaned_df %>% 
    filter(between(length, left = 20000, right = 70000)) %>% 
    select(accession) %>% 
    filter(!str_detect())

fwrite(accessions, "plasmid_accessions_to_download.txt", col.names = F)
    

    
    
