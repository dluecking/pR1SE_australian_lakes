# author: dlueckin
# date: Mon Aug 14 11:58:01 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)
library(rentrez)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# import data -------------------------------------------------------------

big_df <- data.table()

for(file in list.files("hmmsearch_out/", full.names = TRUE, pattern = "domtblout")){
    is_not_empty <- as.numeric(str_split(system(command = paste0("wc -l ", file), intern = TRUE), " ")[[1]][1]) > 13
    
    if(is_not_empty){
        domtbl <- fread(file, skip = 3, fill = TRUE) %>% 
            select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
            filter(V1 != "#")
        
        names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
        
        current_ORF <- str_extract(file, pattern = "ORF[^\\_]*")
        
        # calc alignment length
        domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start
        
        # prepare data for later filtering
        domtbl$query_length <- as.numeric(domtbl$query_length)
        domtbl$score <- as.numeric(domtbl$score)
        
        domtbl <- domtbl %>% 
            filter(score >= 50)
        
        big_df <- rbind(big_df, domtbl)
    }
}



# use function we used earlier to retrieve organism ID --------------------

prot_df <- big_df %>% 
    select(hit, query, score)

names(prot_df) <- c("acc", "hit", "score")

prot_df$organism_acc <- NA
prot_df$organism_name <- NA
prot_df$orf <- str_extract(string = prot_df$hit, pattern = "ORF[^\\.]*")

# Specify your email address for the Entrez API
email <- "dom.luecking@gmail.com"
API_KEY = "ed38a68e4cac02507e4bc585e8913bab5a08"
MAX_ITERATION = 5

# function to retrieve the ID and org name
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


fwrite(prot_df, "prot_df_with_acc_and_id.csv")
