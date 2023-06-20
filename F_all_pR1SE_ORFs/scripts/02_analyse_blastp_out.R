# author: dlueckin
# date: Fri Jun  2 10:43:59 2023

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


# import data -------------------------------------------------------------

blast_df <- fread("blastp_pR1SE_vs_all_complete_genes.out")
blast_df$ORF <- as.numeric(str_extract(blast_df$V1, pattern = "\\d*$"))
count_df <- data.table(table(blast_df$ORF))
count_df$V1 <- as.numeric(count_df$V1)
count_df <- count_df %>% 
    arrange(V1)

tot_contigs <- length(list.files("../../E_check_synteny/genes/"))

ggplot(count_df, aes(x = V1, y = N)) +
    geom_bar(stat = 'identity', color = "black", fill = "grey", alpha = 0.8) +
    theme_minimal() +
    ylim(c(0, tot_contigs)) +
    xlab("ORF#") +
    ylab(paste0("Occurances (total = ", tot_contigs, ")"))

ggsave(plot = last_plot(), filename = "../../plots/all_ORF_occurances.pdf", height = 4, width = 7)
