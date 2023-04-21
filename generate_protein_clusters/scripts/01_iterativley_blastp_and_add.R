# author: dlueckin
# date: Tue Apr 18 14:43:29 2023

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


# functions ---------------------------------------------------------------
# function that takes an input file and blasts against NR, returns the blast as outfmt 6


# function that retrieves a list of accessions and retrieves the sequence to that and writes to a given file
getSequences <- function(accessions, output_file){
    sequences <- entrez_fetch(db="protein", id=accessions, rettype = "fasta")
    write(x = sequences, file = output_file)
}

# iterative search --------------------------------------------------------

ORF <- "ORF_6"
THREADS <- 16

for(ORF in ORFs){
    
    
    we_still_have_a_blast_result = TRUE
    iteration_counter = 0
    
    while(we_still_have_a_blast_result){
        if(iteration_counter == 0){
            input_proteins <- paste0("../initial_proteins/", ORF, ".faa")
        }else{
            input_proteins <- paste0("iterative_blastp/", ORF, "/", ORF, "it_", iteration_counter, "/.faa")
        }
        
        # blast initial proteins
        blast_cmd <- paste0("blastp -query ", input_proteins, " -db ~/bioinf/dbs/nr/nr -outfmt 6 -num_threads ", THREADS)
        blast_result <- system(blast_cmd, intern = TRUE)
        # did we find something new?
        
        # filter the blast df
        
        # download and add proteins to file
        
        # 
        
        
    }
    
}











library(seqinr)
library(rBLAST)

file <- "../initial_proteins/ORF10_initial_proteins.faa"













