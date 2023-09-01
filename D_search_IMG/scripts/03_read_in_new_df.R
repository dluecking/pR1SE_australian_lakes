# author: dlueckin
# date: Wed Aug 30 16:24:30 2023

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



# read df -----------------------------------------------------------------

big_df <- fread("identified_viruses_with_info.tsv")
