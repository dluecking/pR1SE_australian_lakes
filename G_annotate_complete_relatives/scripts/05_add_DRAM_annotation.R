# author: dlueckin
# date: Tue Jun 27 15:14:33 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(plyr) # for rbind.fill

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# read dfs ----------------------------------------------------------------

prot_df <- fread("improved_gene_df.tsv")

# read DRAM  output -------------------------------------------------------

DRAM_df <- data.table()

for(dir in list.dirs("../DRAM_out", recursive = FALSE)){
    tmp_df <- fread(paste0(dir, "/annotations.tsv")) %>% 
        filter(rank %in% c("A", "B", "C"))
        
    DRAM_df <- rbind.fill(DRAM_df, tmp_df)
}

# fix id
DRAM_df$id <- str_remove(DRAM_df$V1, pattern = paste0("^", DRAM_df$fasta, "\\_"))

# join dfs
prot_df$DRAM_kegg <- DRAM_df$kegg_hit[match(prot_df$id, DRAM_df$id)]
prot_df$DRAM_pdfam <- DRAM_df$pfam_hits[match(prot_df$id, DRAM_df$id)]

# write
fwrite(prot_df, file = "improved_gene_df_with_DRAM.tsv")
