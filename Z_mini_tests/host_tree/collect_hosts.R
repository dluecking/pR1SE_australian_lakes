# author: dlueckin
# date: Thu Sep 21 10:31:15 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)

# -------------------------------------------------------------------------

all_files <- list.files("../non_core_region_origin/all_halobacteria/data", full.names = TRUE)
host_accessions <- fread("list_of_host_accessions.txt", header = F)

# collect the hosts -------------------------------------------------------

for(acc in host_accessions$V1){
    file <- all_files[str_detect(all_files, pattern = acc)]
    cmd <- paste0("cp ", file, " host_genomes/")
    system(command = cmd)
}
