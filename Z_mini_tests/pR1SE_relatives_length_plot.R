# author: dlueckin
# date: Wed Jul 26 14:35:54 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)
library(cowplot)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))



# plot length of circular plasmids ----------------------------------------

region_df %>% filter(genome_segment == "circular", integrated_state != "main_chromosome") %>%
    ggplot(aes(x = reorder(contig, sequence_length), y = sequence_length)) +
    geom_bar(stat = 'identity', color = "black", fill = "grey", alpha = 0.8) +
    xlab("Contig") +
    ylab("Sequence Length (bp)") +
    coord_flip() +
    theme_cowplot()


ggsave(plot = last_plot(), file = "../plots/circular_plasmid_length.png", height = 5, width = 6)
