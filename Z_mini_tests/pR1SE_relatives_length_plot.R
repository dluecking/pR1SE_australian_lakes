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



# plot length of region ---------------------------------------------------

region_df %>% filter(genome_segment == "circular", integrated_state != "main_chromosome") %>%
    ggplot(aes(x = reorder(contig, sequence_length), y = region_length_bp)) +
    geom_bar(stat = 'identity', color = "black", fill = "grey", alpha = 0.8) +
    xlab("Contig") +
    ylab("Sequence Length (bp)") +
    coord_flip() +
    theme_cowplot()


ggsave(plot = last_plot(), file = "../plots/circular_plasmid_region_length.png", height = 5, width = 6)



# plot both things in one plot --------------------------------------------

region_df$length_minus_region <- region_df$sequence_length - region_df$region_length_bp

to_plot <- region_df %>% 
    filter(genome_segment == "circular", integrated_state != "main_chromosome") %>% 
    select(contig, length_minus_region, region_length_bp, sequence_length)



a <- gather(to_plot, key = "which", "length", -contig)

# Ensure that 'which' is a factor variable
a$which <- as.factor(a$which)

# Create a new column for ordering by region_length_bp
a <- a %>%
    group_by(contig) %>%
    mutate(order = ifelse(which == "sequence_length", -length, length)) %>%
    ungroup()

# Order the contig variable
a$contig <- factor(a$contig, levels = unique(a$contig[order(a$order)]))

a <- a %>% filter(which != "sequence_length")

ggplot(a, aes(x = contig, y = length, fill = which)) +
    geom_bar(stat = 'identity', color = "black", alpha = 0.8) +
    coord_flip() +
    xlab("Contig") +
    ylab("Sequence Length (bp)") +
    labs(fill = "") +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("length_minus_region" = "grey", "region_length_bp" = "orange"),
                      labels = c("full pR1SE contig", "core region"))

ggsave(plot = last_plot(), file = "../plots/ciruclar_plasmids_region_and_full_contig_dualplot.png", height = 5, width = 6)




# plot region ratio -------------------------------------------------------


# quick calculation
df <- region_df %>% 
    filter(genome_segment == "circular", integrated_state != "main_chromosome") %>% 
    select(contig, length_minus_region, region_length_bp, sequence_length)
df$ratio <- df$sequence_length / df$region_length_bp



ggplot(df, aes(x = reorder(contig, sequence_length), y = ratio)) +
    geom_bar(stat = 'identity', color = "black", fill = "grey", alpha = 0.8) +
    xlab("Contig") +
    ylab("Sequence Length / Core region length ratio") +
    coord_flip() +
    theme_cowplot()
