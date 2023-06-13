# author: dlueckin
# date: Wed May 10 11:50:29 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(rentrez)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# Specify your email address for the Entrez API
email <- "dom.luecking@gmail.com"


# function ----------------------------------------------------------------

getOrganism <- function(acc){
    
    # Use tryCatch to handle errors
    tryCatch({
        
        # Use entrez_summary() to retrieve the summary for the protein accession
        summary_result <- entrez_summary(db="protein", id=acc)
        # Sys.sleep(0.2)
        # Extract the organism name from the summary_result
        org_name <- summary_result["organism"]$organism
        
        # return the organism name
        return(org_name)
        
    }, error = function(e) {
        
        # If an error occurs, return "NA"
        return(NA)
        
    })
    
}

# load ORF clusters -------------------------------------------------------

prot_df <- data.table()

for(file in list.files("../../A_generate_protein_clusters/FINAL_proteins/", full.names = TRUE)){
    accs <- unlist(getName(read.fasta(file = file)))
    accs <- accs[!accs == "Consensus"]
    
    organisms <- sapply(accs, getOrganism)
    
    tmp_df <- data.table(acc = accs,
                         organism = organisms,
                         orf = str_extract(string = file, pattern = "ORF[^\\.]*"))
    prot_df <- rbind(prot_df, tmp_df)
}


fwrite(prot_df, file = "prot_df.tsv")

