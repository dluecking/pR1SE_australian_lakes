# author: dlueckin
# date: Thu Sep 21 14:27:15 2023

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




# load tree ---------------------------------------------------------------

tree <- read.newick(file = "consensus_iq_tree_in_newick.nwk", node.label = "support")
tree <- as.phylo(tree)


# load reigon df and connect the same host info as we did with the --------

# we need family information
region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))

feature_table <- data.table(label = tree$tip.label, host = "")
# ugly regex to be able to match
feature_table$label <- str_remove(sub("^(.*?_.*?_).*", "\\1", feature_table$label), "\\_$")
feature_table$host <- region_df$family[match(feature_table$label, region_df$host_genome_accession)]

# Define a color palette mapping each unique host value to a color
color_palette <- c("#ffbe0b", "#fb5607", "#ff006e", "#8338ec", "#3a86ff", "grey", "darkgreen")
# Add the "color" column to the feature_table based on the host value
feature_table$color[feature_table$host == "Haloarculaceae"] <- "#3a86ff"
feature_table$color[feature_table$host == "Halobacteriaceae"] <- "#fb5607"
feature_table$color[feature_table$host == "Halorubraceae"] <- "#8338ec"
feature_table$color[feature_table$host == "Natrialbaceae"] <- "#ff006e"
feature_table$color[feature_table$host == "Natronoarchaeaceae"] <- "#ffbe0b"



# change tip label to short genome acc
tree$tip.label <- str_remove(sub("^(.*?_.*?_).*", "\\1", tree$tip.label), "\\_$")
# replace label with MGE contig name
tree$tip.label <- region_df$contig[match(tree$tip.label, region_df$host_genome_accession)]


ggtree(tree, branch.length = "none") + 
    geom_tiplab(color = "black", offset = -10) +
    geom_tippoint(color = "black", shape = 23, size = 4, fill = feature_table$color) +
    theme(
        panel.grid = element_line(color = "transparent"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ) +
    scale_x_reverse()


ggsave(plot = last_plot(), file = "../../plots/host_tree.png", height = 9, width = 5, bg = "transparent")
ggsave(plot = last_plot(), file = "../../plots/host_tree.svg", height = 9, width = 5, bg = "transparent")
