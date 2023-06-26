# author: dlueckin
# date: Thu Jun 22 14:15:48 2023

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


# load prot_df ------------------------------------------------------------

prot_df <- fread("prot_df.tsv")

contigs_to_download <- unique(prot_df$organism_acc)
contigs_to_download <- unique(str_remove(contigs_to_download, "^NZ\\_"))



# download fasta ----------------------------------------------------------

for(acc in contigs_to_download){
    seq <- entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
    write(file = paste0("../known_hosts/", acc, ".fasta"), seq)
}
