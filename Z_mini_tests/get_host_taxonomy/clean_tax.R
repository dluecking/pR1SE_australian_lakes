# author: dlueckin
# date: Wed Jul  5 15:02:02 2023

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


# read original df --------------------------------------------------------

tax_df <- fread("tax_df.csv")

tax_df$kingdom <- ""
tax_df$phylum <- ""
tax_df$clade <- ""
tax_df$class <- ""
tax_df$order <- ""
tax_df$family <- ""
tax_df$genus <- ""

for(i in 1:nrow(tax_df)){
    if(tax_df$taxonomy[i] == "")
        next
    
    tax_df$kingdom[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "kingdom_[^,]*"), 
                                    pattern = "kingdom\\_")
    tax_df$phylum[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "phylum_[^,]*"), 
                                   pattern = "phylum\\_")
    tax_df$clade[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "clade_[^,]*"), 
                                  pattern = "clade\\_")
    tax_df$class[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "class_[^,]*"), 
                                  pattern = "class\\_")
    tax_df$order[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "order_[^,]*"), 
                                  pattern = "order\\_")
    tax_df$family[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "family_[^,]*"), 
                                   pattern = "family\\_")
    tax_df$genus[i] <- str_remove(str_extract(tax_df$taxonomy[i], pattern = "genus_[^,]*"), 
                                  pattern = "genus\\_")
}

fwrite(tax_df, "tax_df_cleaned.csv")
