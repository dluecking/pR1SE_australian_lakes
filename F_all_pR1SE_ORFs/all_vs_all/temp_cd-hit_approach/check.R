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

prot_df <- fread("../../../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv")



# read clusters -----------------------------------------------------------

# read original clean df
clean <- fread("combined_relevant_proteins_clustered_15.faa.clstr_CLEAN", header = FALSE)
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
    select(cluster, cluster_annot, contig, presence) %>% 
    distinct()



# remove clusters with less than 4 genes belonging to it ------------------

cluster_df <- as.data.table(table(clean$cluster))
clean$occurances <- cluster_df$N[match(clean$cluster, cluster_df$V1)]

# keep every cluster that either has an annotation or more than 4 members
clean <- clean %>% 
    filter(cluster_annot != "" | occurances >= 4)


# change order ------------------------------------------------------------

clean$ORF_numeric <- as.numeric(str_remove(clean$cluster_annot, "ORF"))
clean <- arrange(clean, ORF_numeric, desc(occurances))

clean$order <- 0
current_cluster <- clean$cluster[1]
current_order <- 1

for(i in 1:nrow(clean)){
    # if we are in the same cluster
    if(clean$cluster[i] == current_cluster)
        clean$order[i] <- current_order
    else{
        current_cluster <- clean$cluster[i]
        current_order <- current_order + 1
        
        clean$order[i] <- current_order
    }

}
    
    

# y order, so we have KX on top, then alphabetically ----------------------

contig_order <- sort(unique(clean$contig))
contig_order <- contig_order[contig_order != "KX906370.1"]
contig_order <- c("KX906370.1", contig_order)
    
# first attempt to plot ---------------------------------------------------

ggplot(clean, aes(cluster, contig, fill = cluster_annot, alpha = presence)) + 
    geom_tile(colour = "black") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          text = element_text(size = 9),
          axis.text.x = element_blank()) +
    scale_x_discrete(limits=(unique(clean$cluster))[order(unique(clean$order))],
                     position = "top") +
    scale_y_discrete(limits=rev(contig_order)) +
    xlab("") +
    ylab("") +
    theme(legend.position = "None") +
    scale_fill_manual(values = c("ORF6" = "#f94144",
                                 "ORF8" = "#f3722c",
                                 "ORF10" = "#f8961e",
                                 "ORF17" = "#f9c74f",
                                 "ORF21" = "#90be6d",
                                 "ORF23" = "#43aa8b",
                                 "ORF24" = "#577590",
                                 "unknown" = "#FBFEF9"))
    
ggsave(plot = last_plot(), filename = "../../../plots/protein_cluster_heatmap.svg",
       height = 7, width = 11)


# for plotting in biorender
u_df <- clean %>% 
    filter(cluster_annot == "")
