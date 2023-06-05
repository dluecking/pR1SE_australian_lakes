# author: dlueckin
# date: Wed May 10 15:18:46 2023

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



# import df ---------------------------------------------------------------

prot_df <- fread("prot_df.tsv")
protein_accession <- "WP_004048516"
taxid <- 2248

get_genome_accessions <- function(taxid) {
    query <- paste0("txid", taxid, "[Organism:exp]")
    esearch <- entrez_search(db = "nucleotide", term = query)
    id_list <- esearch$ids
    esummary <- entrez_summary(db = "nucleotide", id = id_list)
    genbank_accessions <- esummary$caption[grepl("genbank", esummary$caption)]
    genome_accessions <- sub("\\s+.*", "", genbank_accessions)
    return(genome_accessions)
}










#  old --------------------------------------------------------------------


host_df <- fread("known_hosts_df.tsv")

acc <- "WP_004048516"

# get protein links
protein_links <- entrez_link(dbfrom='protein', id = acc, db='all')

# use taxid number to get taxid links
tax_id <- protein_links$links$protein_taxonomy
tax_links <- entrez_link(dbfrom='taxonomy', id = tax_id, db='all')

# get genome based on taxid
genome_id <- tax_links$links$taxonomy_genome

genome_id <- "448447358"
genome_fasta <- entrez_fetch(db = "nuccore", id = genome_id, rettype = "fasta")





download_genome <- function(genome_accession, file_name) {
    # set up NCBI API parameters
    db <- "nuccore"
    id <- "NZ_AOJE01000051.1"
    rettype <- "gb"
    
    # construct the download URL
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=", db, "&id=", id, "&rettype=", rettype)
    
    # download the file and save it
    download.file(url, destfile = "here", mode = "wb")
}
