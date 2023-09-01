# author: dlueckin
# date: Tue Aug  8 17:21:15 2023

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


# load protein df ---------------------------------------------------------

prot_df <- fread("../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv")


# load clean_df -----------------------------------------------------------

clean <- fread("../F_all_pR1SE_ORFs/all_vs_all/temp_cd-hit_approach/combined_relevant_proteins_clustered_15.faa.clstr_CLEAN", header = FALSE)

# clean up the clean_df
clean$V1 <- str_remove(clean$V1, pattern =">")
clean$annot <- prot_df$annot[match(clean$V1, prot_df$id)]

# preserve order in case you fuck up and sort
clean$order <- c(1:nrow(clean))

# add cluster
clean$cluster <- ""

current_cluster <- "NIX"

for(i in 1:nrow(clean)){
    if(startsWith(clean$V1[i], prefix = "---")){
        current_cluster <- clean$V1[i]
        next                      
    }
    clean$cluster[i] <- current_cluster
}

# remove the cluster rows
clean <- clean %>% 
    filter(cluster != "")


# add a cluster wide annotation (this changes most likely not much, but some!)
clean$cluster_annot <- ""
for(i in 1:nrow(clean)){
    tmp_df <- clean %>% 
        filter(cluster == clean$cluster[i]) %>% 
        filter(str_detect(annot, "ORF")) %>% 
        select(annot)
    
    if(nrow(tmp_df) == 0){
        clean$cluster_annot[i] <- ""
    }else{
        clean$cluster_annot[i] <- tmp_df %>% 
            slice(1) %>% 
            pull()
    }
}
rm(i, tmp_df, current_cluster)



# some cosmetics ----------------------------------------------------------
# add presence column
clean$presence <- 1

# add contig
clean$contig <- str_remove(clean$V1, pattern = "\\_\\d*$")
clean$contig <- str_remove(clean$contig, pattern = "\\_length.*$")

# remove ---
clean$cluster <- str_remove(clean$cluster, pattern = "--- ")

# remove .faa
clean$cluster_annot <- str_remove(clean$cluster_annot, pattern = "\\.faa$")

# select the relevant contigs AND
# if a contig has two hits in the same cluster, we count only one -> distinct
clean <- clean %>% 
    select(cluster, cluster_annot, contig, presence, V1) %>% 
    distinct()

names(clean) <- c("cluster", "cluster_annotation", "contig", "presence", "gene_id")


# remove clusters with less than 4 genes belonging to it ------------------

cluster_df <- as.data.table(table(clean$cluster))
clean$occurances <- cluster_df$N[match(clean$cluster, cluster_df$V1)]


# save df -----------------------------------------------------------------

fwrite(clean, "clean_cluster_df.csv")







# quick: figure out what cluster 14 and cluster 16 are --------------------

c14 <- clean %>% 
    filter(cluster == "Cluster 14")
c14$nr_annot <- prot_df$nr_annotation[match(c14$gene_id, prot_df$id)]

fwrite(c14, "cluster_14_annot.csv")

c16 <- clean %>% 
    filter(cluster == "Cluster 16")
c16$nr_annot <- prot_df$nr_annotation[match(c16$gene_id, prot_df$id)] 

fwrite(c16, "cluster_16_annot.csv")









