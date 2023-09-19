# author: dlueckin
# date: Thu Aug 10 10:20:42 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(vegan)
# library(hclust)
library(ggplotify)
library(ggtree)
library(ape)
library(googlesheets4)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# import data -------------------------------------------------------------

cluster_df <- fread("clean_cluster_df.csv") %>% 
    filter(cluster_annotation != "" | occurances >= 4)

cluster_df$final_cluster <- cluster_df$cluster
cluster_df$final_cluster[cluster_df$cluster_annotation != ""] <- cluster_df$cluster_annotation[cluster_df$cluster_annotation != ""]


# create new wide df ------------------------------------------------------

wide <- data.frame(cluster = unique(cluster_df$final_cluster))
wide[unique(cluster_df$contig)] <- 0

# fill it
for(row in 1:nrow(wide)){
    for(column in 2:length(wide)){
        value <- cluster_df %>% 
            filter(final_cluster == wide$cluster[row]) %>% 
            filter(contig == names(wide)[column]) %>% 
            nrow()
        
        if(value >= 1)
            wide[row,column] <- 1 
        
    }
}
rm(row, column)



# prepare for tree --------------------------------------------------------

wide_transposed <- t(wide %>% select(-cluster))
dist <- vegdist(wide_transposed)

hclust_avg <- hclust(dist, method = "average")
p_classic <- plot(hclust_avg, xlab = "", ylab = "distance",
                  hang = -1, direction = "left")
# not great.  I have to export it. Or figure out





# we need family information
region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
region_df <- region_df %>% 
    select(contig, family)

feature_table <- data.table(label = tree$tip.label, host = "")
feature_table$host <- region_df$family[match(feature_table$label, region_df$contig)]

# Define a color palette mapping each unique host value to a color
color_palette <- c("#ffbe0b", "#fb5607", "#ff006e", "#8338ec", "#3a86ff", "grey", "darkgreen")

# Add the "color" column to the feature_table based on the host value
feature_table[, color := color_palette[match(host, unique(host))]]


tree <- as.phylo(hclust_avg)
plot(tree, type = "tidy", direction = "right", 
     show.node.label = TRUE, 
     no.margin = TRUE, font = 1,
     col = feature_table$color)

# 