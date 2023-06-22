# author: dlueckin
# date: Wed Jun 21 18:07:05 2023

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



# import hmm results ------------------------------------------------------

hmm_df <- fread("01_hmm_results/combined_res.out", header = F, sep =" ")
