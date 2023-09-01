# author: dlueckin
# date: Tue Jul  4 15:39:06 2023

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

prot_df <- fread("improved_gene_df_with_DRAM.tsv")


# import blast out --------------------------------------------------------

blast_df <- fread("../diamond_blastp_nr_out/all_pR1SE_relatives_genes_vs_nr.out")
names(blast_df) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle", " "))


# filter ------------------------------------------------------------------

# remove hypos
blast_df <- blast_df %>% 
    filter(!str_detect(stitle, "hypothetical")) 

blast_df$stitle <- str_remove(blast_df$stitle, blast_df$sseqid)



# add to prot_df ----------------------------------------------------------

prot_df$nr_annotation <- NA
prot_df$nr_bitscore <- NA

for(i in 1:nrow(prot_df)){
    tmp_df <- blast_df %>% 
        filter(qseqid == prot_df$id[i]) %>% 
        slice_max(order_by = bitscore, n = 1)
    
    if(nrow(tmp_df) >= 1){
        prot_df$nr_annotation[i] <- tmp_df$stitle[1] %>% 
            trimws()
        prot_df$nr_bitscore[i] <- tmp_df$bitscore[1]
    }    
}

fwrite(prot_df, "improved_gene_df_with_DRAM_and_nr.tsv")




# quick -------------------------------------------------------------------


prot_df <- fread("improved_gene_df_with_DRAM_and_nr.tsv") %>% filter(annot == "ORF6.faa")
cp_df <- data.table()

for(i in 1:nrow(region_df)){
    tmp_df <- prot_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(`orf#` >= region_df$manual_min_orf[i]) %>% 
        filter(`orf#` <= region_df$manual_max_orf[i])
    
    cp_df <- rbind(cp_df, tmp_df)
}
fwrite(prot_df, "ORF6_df_full.csv")
fwrite(cp_df, "ORF6_df_only_complete.csv")


safe_df <- cp_df %>% 
    filter(aln_length >= 1000)


# get very safe proteins
for()