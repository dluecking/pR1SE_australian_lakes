# author: dlueckin
# date: Fri Sep 22 17:12:01 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(dplyr)


# # working directory -------------------------------------------------------
# this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(this_dir)
# print(paste0("Setting wd to: \n ", this_dir))
# 
# # region_df
# region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
#                                       sheet = "pR1SE_relatives"))

# I've retrieved nuccore2 based on this query:
# "archaeal"[Organism] OR "Archaea"[Organism]) AND "plasmid"[Filter] NOT "shotgun"[All Fields] NOT "MAG"[All Fields] 

# explore -----------------------------------------------------------------

# NCBI search term:
# "archaeal"[Organism] OR "Archaea"[Organism]) AND "plasmid"[Filter] NOT "shotgun"[All Fields] 

# Read the file into a character vector
file_lines <- readLines("nuccore_result2.txt")

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


# clean up df -------------------------------------------------------------
# remove duplicates
plasmid_df <- plasmid_df %>% 
    mutate(accession = str_remove(accession, "NZ\\_")) %>% 
    distinct()


# add class based on accession --------------------------------------------

tax_df <- data.table(id = plasmid_df$accession,
                     tax_id = "",
                     taxonomy = "")

for(i in 1:nrow(tax_df)){
    cmd1 <- paste0("esearch -db nucleotide -query '", tax_df$id[i], "'|esummary|xtract -pattern TaxId -element TaxId")
    t_id <- system(command = cmd1, intern = TRUE)
    if(length(t_id) >= 1){
        tax_df$tax_id[i] <- t_id
    }else{
        tax_df$tax_id[i] <- NA
    }
    if(i %% 10 == 0){
        print(i)
    }
    
}
print("done getting taxId")

for(i in 1:nrow(tax_df)){
    if(!is.na(tax_df$tax_id[i])){
        if(i %% 10 == 0){
            print(i)
        }
        
        cmd2 <- paste0("efetch -db taxonomy -id ", tax_df$tax_id[i], " -format xml | xtract -pattern Taxon -first TaxId -element Taxon -block '*/Taxon' -unless Rank -equals 'no rank' -tab ',' -sep '_' -element Rank,ScientificName")
        tax_df$taxonomy[i] <- system(command = cmd2, intern = TRUE)
    }
}
print("done getting taxonomy")



# combine -----------------------------------------------------------------

plasmid_df$taxonomy <- tax_df$taxonomy[match(plasmid_df$accession, tax_df$id)]
glimpse(plasmid_df)
fwrite(plasmid_df, "plasmid_df_with_tax.csv")