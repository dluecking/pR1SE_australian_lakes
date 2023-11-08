# author: dlueckin
# date: Mon Oct 16 13:25:11 2023

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

# positions ---------------------------------------------------------------

positions <- c(
    "Cell cycle control, cell division, chromosome partitioning",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Posttranslational modification, protein turnover, chaperones",
    "Signal transduction mechanisms",
    "Defense mechanisms",
    "Extracellular structures",
    "Nuclear structure",
    "Cytoskeleton",
    
    "Mobilome: prophages, transposons",
    "Intracellular trafficking, secretion, and vesicular transport",
    
    "Energy production and conversion",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Carbohydrate transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism",
    
    "RNA processing and modification",
    "Chromatin structure and dynamics",
    "Translation, ribosomal structure and biogenesis",
    "Transcription",
    "Replication, recombination and repair",
    
    "General function prediction only",
    "Function unknown"
)

# import data -------------------------------------------------------------

data <- fread("df_for_plotting.csv") 

COG_df <- fread("../../../template_classifier_count.tsv")
COG_df$delta <- data$percentage[match(COG_df$LETTER, data$LETTER) & data$label == "pR1SE"]  -
    data$percentage[match(COG_df$LETTER, data$LETTER) & data$label == "host"]

COG_df$positive <- COG_df$delta >= 0

COG_df$SHORT_DESCRIPTION <- factor(COG_df$SHORT_DESCRIPTION, levels = COG_df$SHORT_DESCRIPTION[order(COG_df$delta)])


ggplot(COG_df, aes(x = SHORT_DESCRIPTION, y = delta, fill = positive)) +
    geom_bar(stat = 'identity', color = "black") +
    theme_cowplot(font_size = 8,) +
    xlab("COG Category") +
    ylab("") +
    coord_flip() +
    scale_fill_manual(values = c("grey", "orange")) +
    theme(legend.position = "None")


ggsave(plot = last_plot(), "../plots/COG_outside_region_vs_hosts.png", height = 3, width = 4)
ggsave(plot = last_plot(), "../plots/COG_outside_region_vs_hosts.svg", height = 3, width = 4)
