# author: dlueckin
# date: Wed May 31 15:23:12 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gggenes)
library(ggridges)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# import df and filter  --------------------------------------------------

gene_df <- fread("cleaned_gene_df_complete.tsv")



# plot sizes --------------------------------------------------------------

box <- ggplot(gene_df, aes(x = annot, y = length, fill = annot)) +
    geom_boxplot() +
    theme_minimal() +
    scale_fill_manual(values = c("ORF6" = "#f94144",
                                 "ORF8" = "#f3722c",
                                 "ORF10" = "#f8961e",
                                 "ORF17" = "#f9c74f",
                                 "ORF21" = "#90be6d",
                                 "ORF23" = "#43aa8b",
                                 "ORF24" = "#577590",
                                 "unknown" = "#FBFEF9")) +
    xlab("") +
    ylab("Length (bp)") +
    theme(legend.position = "None") +
    ggtitle("Avg ORF size (bp)")
ggsave(plot = box, "../../plots/average_ORF_lengths_boxplot.pdf", height = 4, width = 7)



box2 <- ggplot(gene_df, aes(x = annot, y = score, fill = annot)) +
    geom_boxplot() +
    theme_minimal() +
    scale_fill_manual(values = c("ORF6" = "#f94144",
                                 "ORF8" = "#f3722c",
                                 "ORF10" = "#f8961e",
                                 "ORF17" = "#f9c74f",
                                 "ORF21" = "#90be6d",
                                 "ORF23" = "#43aa8b",
                                 "ORF24" = "#577590",
                                 "unknown" = "#FBFEF9")) +
    xlab("") +
    ylab("Score") +
    theme(legend.position = "None") +
    ggtitle("Avg ORF hmm bitscore")
ggsave(plot = box, "../../plots/average_ORF_score_boxplot.pdf", height = 4, width = 7)
