# author: dlueckin
# date: Wed May 31 13:55:57 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gggenes)



# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# import gene data --------------------------------------------------------

gene_df <- data.table()

for(file in list.files("../sco/complete/")){
    tmp_df <- fread(paste0("../sco/complete/", file), skip = 2, sep = "_")
    tmp_df$contig <- str_remove(file, "\\.fasta\\.sco")
    names(tmp_df) <- c("orf#", "start", "end", "strand", "contig")
    
    gene_df <- rbind(gene_df, tmp_df)
}
rm(tmp_df, file)

gene_df$strand[gene_df$strand == "+"] <- TRUE
gene_df$strand[gene_df$strand == "-"] <- FALSE
gene_df$`orf#` <- as.numeric(str_remove(gene_df$`orf#`, ">"))


# import hmm results ------------------------------------------------------

hmm_df <- data.table()

for(file in list.files("../hmm_results/complete", full.names = TRUE, pattern = "domtblout.txt")){
    is_not_empty <- as.numeric(str_split(system(command = paste0("wc -l ", file), intern = TRUE), " ")[[1]][1]) > 13
    
    if(is_not_empty){
        domtbl <- fread(file, skip = 3, fill = TRUE, header = FALSE) %>% 
            select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
            filter(V1 != "#")
        
        names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")

        # calc alignment length
        domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start
        
        # prepare data for later filtering
        domtbl$query_length <- as.numeric(domtbl$query_length)
        domtbl$score <- as.numeric(domtbl$score)
        
        hmm_df <- rbind(hmm_df, domtbl)
    }
}
rm(domtbl, is_not_empty, file)


# combining the plots -----------------------------------------------------

gene_df$id <- paste0(gene_df$contig, "_", str_remove(gene_df$`orf#`, ">"))
hmm_df$id <- str_remove(hmm_df$hit, ">")

gene_df$annot <- hmm_df$query[match(gene_df$id, hmm_df$id)]
gene_df$score <- hmm_df$score[match(gene_df$id, hmm_df$id)]
gene_df$evalue <- hmm_df$evalue[match(gene_df$id, hmm_df$id)]
gene_df$aln_length <- hmm_df$aln_length[match(gene_df$id, hmm_df$id)]




# filter data -------------------------------------------------------------

ORF_WINDOW_SIZE = 5
gene_df_filtered <- data.table()

for(acc in unique(gene_df$contig)){
    tmp_df <- gene_df %>% 
        filter(acc == contig) %>% 
        filter(annot != "") 
    
    # this gives us the ORF window
    
    maximum_orf_number_overall <- gene_df %>% 
        filter(acc == contig) %>% 
        select(`orf#`)
    maximum_orf_number_overall <- max(maximum_orf_number_overall$`orf#`)
    
    minimum_orf <- max(0, min(tmp_df$`orf#`)-ORF_WINDOW_SIZE)
    maximum_orf <- min(maximum_orf_number_overall, max(tmp_df$`orf#` + 5))
    
    # now filter the original df based on this
    t_df <- gene_df %>% 
        filter(contig == acc) %>% 
        filter(`orf#` >= minimum_orf, `orf#` <= maximum_orf)
    gene_df_filtered <- rbind(gene_df_filtered, t_df)
}
rm(t_df, acc, tmp_df, maximum_orf_number_overall)


# this is correct, but there is plenty of space between the orfs in most cases...especially "NZ_CP101161.1"
fwrite(gene_df_filtered, "gene_df.tsv")



