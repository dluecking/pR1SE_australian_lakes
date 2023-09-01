# author: dlueckin
# date: Wed Jul  5 11:13:34 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(dplyr)

# working directory -------------------------------------------------------
# this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(this_dir)
# print(paste0("Setting wd to: \n ", this_dir))

pR1SE_df <- fread("all_contigs.csv")

tax_df <- data.table(id = pR1SE_df$contig,
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
    
}

for(i in 1:nrow(tax_df)){
    if(!is.na(tax_df$tax_id[i])){
        print(i)
        print(tax_df$tax_id[i])
        
        cmd2 <- paste0("efetch -db taxonomy -id ", tax_df$tax_id[i], " -format xml | xtract -pattern Taxon -first TaxId -element Taxon -block '*/Taxon' -unless Rank -equals 'no rank' -tab ',' -sep '_' -element Rank,ScientificName")
        tax_df$taxonomy[i] <- system(command = cmd2, intern = TRUE)
    }
}
fwrite(tax_df, "tax_df.csv")