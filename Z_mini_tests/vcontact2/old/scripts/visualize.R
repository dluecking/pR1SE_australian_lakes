# author: dlueckin
# date: Wed Sep 13 16:59:31 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)
library(igraph)
library(ggraph)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
region_df <- region_df %>% 
    filter(contig_category == "complete")


# load network file --------------------------------------------------------

network_df <- fread("vcontact2_out/c1.ntw")

network_df <- network_df %>% 
    sample_n(1000)


network_df$is_pR1SE <- if_else(network_df$V1 %in% region_df$contig, "YES", "NO")

network_data <- graph_from_data_frame(network_df, directed = F)


ggraph(network_data, layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(color = network_df$is_pR1SE)) +
    geom_node_text(aes(label = "")) +
    theme_void()
