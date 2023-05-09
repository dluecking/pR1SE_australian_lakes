# author: dlueckin
# date: Wed May  3 10:00:01 2023

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


# import ------------------------------------------------------------------
data <- fread("pR1SE core set - interpro_out.tsv")

files <- list.files("../05_curated_blastp_proteins/", full.names = TRUE)
for(i in 1:length(files)){
    accs <- getName(read.fasta(files[i]))
    for(acc in accs){
        orf <- str_extract(string = files[i], pattern = "ORF[^\\_]*")
        data$ORF[data$V1 == acc] <- orf
            
    }
}
data <- setorder(data,ORF)
fwrite(data, file = "sorted_interpro_out.tsv")
