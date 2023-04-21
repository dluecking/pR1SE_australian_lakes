# author: dlueckin
# date: Tue Apr 18 10:54:15 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(tidyr)
library(rentrez)


# load clean df -----------------------------------------------------------
# local
INPUT_PATH <- "input/pR1SE core set - clean input.csv"

# read input df
input_df <- fread(INPUT_PATH)

# convert to long
input_df <- input_df %>% 
    pivot_longer(-c(pR1SE, `predicted function`), names_to = "species", values_to = "accession") %>% 
    mutate(species = gsub(" ", "_", species))
names(input_df) <- c("ORF", "predicted_function", "species", "accession")


# iterate over given ORFs -------------------------------------------------

for(current_ORF in input_df$ORF){
    list_of_accessions <- input_df %>% 
        filter(current_ORF == ORF) %>% 
        select(accession) %>% 
        unlist()
    
    FILENAME <- paste0("initial_proteins/", current_ORF, "_initial_proteins.faa")
    
    sequences <- entrez_fetch(db="protein", id=list_of_accessions, rettype = "fasta")
    write(x = sequences, file = FILENAME)
}






