# author: dlueckin
# date: Mon Jun 26 16:58:19 2023

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

out <- fread("interpro_out.tsv", header = F, sep = "\t")
out <- out %>% 
    filter(V9 != "-")
out <- out %>% 
    filter(as.numeric(V9) <= 0.00001)
out$orf <- as.numeric(str_extract(out$V1, pattern = "\\d*$"))
out <- arrange(out, orf)

fwrite(out, "KX687704_interpro_out_clean.csv", quote = TRUE)
