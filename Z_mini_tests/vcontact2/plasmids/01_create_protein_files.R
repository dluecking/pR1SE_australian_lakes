# author: dlueckin
# date: Wed Jun 28 11:02:58 2023

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


# genes for the pR1SE were just copied from the "virus only" approach in the parent folder
# (vs2) [11:11:16] ~/.../vcontact2/plasmids/genes $ cp ../../input_files/protein_file.faa .



# the plasmids gene file was created by prodigal:
# (vs2) [11:07:58] ~/.../vcontact2/plasmids/genes $ prodigal -i archaeal_plasmids_between_20k_70k.fasta -a ../genes/archaeal_plasmids_between_20k_70k_genes.faa -p meta

# files were combined
# (vs2) [11:12:35] ~/.../vcontact2/plasmids/genes $ cat * > combined_protein_file.faa
