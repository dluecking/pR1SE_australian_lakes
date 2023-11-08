# author: dlueckin
# date: Wed Aug  9 17:27:14 2023

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



# load clean-df -----------------------------------------------------------

clean_df <- fread("clean_cluster_df.csv")



# filter less than 4 hits -------------------------------------------------

clean_df <- clean_df %>% 
    filter(cluster_annotation != "" | occurances >= 4)



# collapse the known ORFs into 1 cluster ----------------------------------

clean_df$final_cluster <- clean_df$cluster
clean_df$final_cluster[clean_df$cluster_annotation != ""] <- clean_df$cluster_annotation[clean_df$cluster_annotation != ""]



# set order ---------------------------------------------------------------
# determine the order of the x axis: first the ORFs, then the cluster in descending order
x_order <- c("ORF6", "ORF8", "ORF10", "ORF17", "ORF21", "ORF23", "ORF24" )
x_order <- c(x_order, clean_df %>% 
                 filter(str_detect(final_cluster, "Cluster")) %>% 
                 arrange(desc(occurances)) %>% 
                 select(final_cluster) %>% 
                 unique() 
             %>% pull())

# contig_order for now is KX on top, then the rest, but will change with tree!
contig_order <- sort(unique(clean_df$contig))
contig_order <- contig_order[contig_order != "KX687704.1"]
contig_order <- c("KX687704.1", contig_order)

ggplot(clean_df, aes(final_cluster, contig, fill = cluster_annotation, alpha = presence)) + 
    geom_tile(colour = "black") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_blank()) +
    scale_x_discrete(limits=x_order, position = "top") +
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

ggsave(last_plot(), filename = "../plots/FIGURE_1_pangenome_collapsed.svg", height = 6, width = 16)
ggsave(last_plot(), filename = "../plots/FIGURE_1_pangenome_collapsed.png", height = 6, width = 16)
