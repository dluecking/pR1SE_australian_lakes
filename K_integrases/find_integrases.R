# author: dlueckin
# date: Thu Jun 29 12:41:08 2023

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


search_terms <- paste(c("integrase", "recombinase", "intergrase"), collapse = "|")
search_region <- 20

# read google sheet -------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df$integrase_orf <- ""
manual_region_df$integrase_annotation <- ""

# read big interpro_out df-------------------------------------------------

prot_df <- fread("../G_annotate_complete_relatives/scripts/interpro_out_full.tsv")
prot_df$annotation <- tolower(prot_df$annotation)


# add integrase position and annotation to df -----------------------------

for(i in 1:nrow(manual_region_df)) {
    ORF_POSITIONS <- vector()
    INTEGRASE_ANNOTATION <- vector()
    
    if(is.na(manual_region_df$manual_region_length[i])){
        next
    }
    
    acc <- manual_region_df$contig[i]
    int_df <- prot_df %>% 
        filter(contig == acc) %>% 
        filter(str_detect(pattern = search_terms,
                          string = annotation))
    
    # if not empty
    if(nrow(int_df) >= 1){
        # check ech integrase if its close enough
        for(j in 1:nrow(int_df)){
            if(between(int_df$orf[j], 
                       left = manual_region_df$manual_min_orf[i] - search_region,
                       right = manual_region_df$manual_max_orf[i] + search_region)){
                
                # check if the ORF already is in the list
                # this could be the case, if there are multiple annotations for each integrase
                if(!int_df$orf[j] %in% ORF_POSITIONS){
                    ORF_POSITIONS <- append(ORF_POSITIONS, int_df$orf[j])
                    # if we have multiple, get all annotations
                    best_annotation <- int_df %>% 
                        filter(int_df$orf == int_df$orf[j]) %>%
                        select(annotation) %>% 
                        pull()
                    best_annotation <- paste(best_annotation, collapse = "--")
                    INTEGRASE_ANNOTATION <- append(INTEGRASE_ANNOTATION, best_annotation)
                    
                }
                    
            }
        }
    }
    # fill in
    manual_region_df$integrase_orf[i] <- paste(ORF_POSITIONS, collapse = ",")
    manual_region_df$integrase_annotation[i] <- paste(INTEGRASE_ANNOTATION, collapse = ",")
}

fwrite(manual_region_df %>% select(contig, integrase_orf, integrase_annotation), file = "integrase_positions_and_annotations.csv")


# get integrase AA sequence and write to file -----------------------------

for(i in 1:nrow(manual_region_df)){
    # only if we have some
    if(unlist(manual_region_df$integrase_orf[i]) == "")
        next
    
    acc <- manual_region_df$contig[i]
    integrase_positions <- unlist(str_split(unlist(manual_region_df$integrase_orf[i]), ","))
    integrase_accessions <- unlist(paste0(acc, "_", integrase_positions))
    
    # read gene file
    gene_file <- paste0("../E_check_synteny/genes/", acc, ".fasta.faa")
    all_seqs <- read.fasta(gene_file)
    relevant_seqs <- all_seqs[integrase_accessions]
    
    # write to file
    write.fasta(relevant_seqs, getName(relevant_seqs), file.out = paste0("integrases/", acc, "_integrases.faa"))
}
