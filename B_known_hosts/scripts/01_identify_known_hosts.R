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
API_KEY = "ed38a68e4cac02507e4bc585e8913bab5a08"
MAX_ITERATION = 10

# function ----------------------------------------------------------------

getOrganismAndNuccoreAcc <- function(acc){
    # Use tryCatch to handle errors
    tryCatch({
        # use entrez link to find the protein_nuccore link
        link_result <- entrez_link(dbfrom = "protein", id = acc, db = "all", api_key = API_KEY)
        Sys.sleep(0.1)        

        protein_nuc <- link_result$links["protein_nuccore"]
        
        # catcht the weird entrez results
        if(is.null(protein_nuc)){
            protein_nuc <- unlist(link_result[[1]]$links["protein_nuccore"])
        }
        if(is.list(protein_nuc)){
            protein_nuc <- unlist(protein_nuc)[1]
        }
        
        # Use entrez_summary() to retrieve the summary for the protein accession
        summary_result <- entrez_summary(db="nuccore", id = protein_nuc,  api_key = API_KEY)
        Sys.sleep(0.1)
        
        # Extract the organism name from the summary_result
        org_name <- paste0(summary_result$organism, "_", summary_result$strain)
        org_acc <- summary_result$caption
        
        # return the organism name
        return(paste0(org_name, "____", org_acc))
        
    }, error = function(e) {
        # If an error occurs, return "NA"
        return(NA)
    })
    
}

# load ORF clusters -------------------------------------------------------
# set up prot_df
prot_df <- data.table()
for(file in list.files("../../A_generate_protein_clusters/FINAL_proteins", full.names = TRUE)){
    accs <- unlist(getName(read.fasta(file = file)))
    accs <- accs[!accs == "Consensus"]
    
    tmp_df <- data.table(acc = accs,
                         organism_name = NA,
                         organism_acc = NA,
                         orf = str_extract(string = file, pattern = "ORF[^\\.]*"))
    prot_df <- rbind(prot_df, tmp_df)
}


# fill acc und org name 
# do this 10 times, and see how many we get
j = 1
while(j <= MAX_ITERATION){
    for(i in 1:nrow(prot_df)){
        if(is.na(prot_df$organism_acc[i])){
            Name_and_Acc <- getOrganismAndNuccoreAcc(prot_df$acc[i])
            organism_name <- str_remove(Name_and_Acc, "\\_\\_\\_\\_.*$")
            organism_acc <- str_remove(Name_and_Acc, "^.*\\_\\_\\_\\_")
            
            prot_df$organism_name[i] <- organism_name
            prot_df$organism_acc[i] <- organism_acc
            
        }
    }
    print(paste0("iteration: ", j))
    j <- j + 1
}


# remove NZ\\_ in order to avoid the stupid doubles
prot_df$organism_acc <- str_remove(prot_df$organism_acc, "NZ\\_")


fwrite(prot_df, file = "prot_df.tsv")

