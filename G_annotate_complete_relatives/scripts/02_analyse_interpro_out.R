# author: dlueckin
# date: Wed Jun 21 11:53:40 2023

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

# interpro_out <- rbindlist(sapply(list.files("../interpro_out", pattern = , full.names = TRUE), fread))


# read in data ------------------------------------------------------------

interpro_out <- data.table()
for(file in list.files("../interpro_out/", full.names = TRUE)){
    tmp_df <- fread(file)
    interpro_out <- rbind(interpro_out, tmp_df)
}


# some cosmetics ----------------------------------------------------------

interpro_out <- interpro_out %>% 
    select(V1, V3, V4, V6, V7, V8, V9)

names(interpro_out) <- c("gene_id", "gene_length", "DB", "annotation", "hit_start", "hit_end", "evalue")

# add contig id
interpro_out$contig <- str_remove(interpro_out$gene_id, ">")
interpro_out$contig <- str_remove(interpro_out$contig, "\\_\\d*$")

# remove "-" and change to numeric
interpro_out[interpro_out$evalue == "-"]$evalue <- 100 
interpro_out$evalue <- as.numeric(interpro_out$evalue)

# remove ">"
interpro_out$gene_id <- str_remove(interpro_out$gene_id, ">")

# add orf numbers
interpro_out$orf <- str_extract(interpro_out$gene_id, "\\d*$")


# remove non hits or low scoring hits -------------------------------------

interpro_out <- interpro_out %>% 
    filter(evalue <= 0.00001) %>% 
    filter(annotation != "-")
fwrite(interpro_out, "interpro_out_full.tsv", sep = "\t")


# read contig_border_df and remove unnceessary annotations ----------------

interpro_out <- interpro_out %>% group_by(gene_id) %>% slice_min(n = 2, evalue)


# combine this DF with the annotation DF from E_synteny -------------------

gene_df <- fread("../../E_check_synteny/scripts/gene_df.tsv")
gene_df$interpro_annot <- interpro_out$annotation[match(gene_df$id, interpro_out$gene_id)]
gene_df$interpro_evalue <- interpro_out$evalue[match(gene_df$id, interpro_out$gene_id)]


fwrite(gene_df, "improved_gene_df.tsv", sep = "\t")

