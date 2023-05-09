# author: dlueckin
# date: Wed May  3 11:18:19 2023

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



# import all blast results ------------------------------------------------

data <- rbindlist(lapply(list.files("../02_blastp_results/", full.names = TRUE), fread)) %>% distinct()
