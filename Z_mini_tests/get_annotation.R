# author: dlueckin
# date: Wed Jun 28 16:42:41 2023

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


library(taxonomizr)
prepareDatabase(sqlFile = "accessionTaxa.sql")
