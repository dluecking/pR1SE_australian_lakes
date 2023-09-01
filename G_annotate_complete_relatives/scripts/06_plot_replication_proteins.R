# author: dlueckin
# date: Thu Jun 29 15:40:41 2023

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

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives")) %>% 
    filter(contig_category == "complete")

to_plot <- as.data.table(table(manual_region_df$replication_group))
colors <- c("#04a777", "#FB8B24", "#5f758e", "grey", "#bcd2ee")

ggplot(to_plot, aes(x = reorder(V1, -N), y = N, fill = V1)) +
    geom_bar(stat = "identity", alpha= 0.8, color = "black") +
    theme_minimal() +
    xlab("Replication Type") +
    ggtitle("Replication Types of complete pR1SE Relatives",
            subtitle = "N = 41") +
    scale_fill_manual(values=colors) +
    theme(legend.position = "None")

ggsave(plot = last_plot(), filename = "../../plots/replication_types.pdf", height = 5, width = 5)
