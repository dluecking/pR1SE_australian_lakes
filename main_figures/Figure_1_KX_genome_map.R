# author: dlueckin
# date: Wed Aug  9 16:46:29 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gggenes)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# read gene df ------------------------------------------------------------

gene_df <- fread("../E_check_synteny/sco/KX687704.1.fasta.sco", skip = 2, sep = "_")
names(gene_df) <- c("orf#", "start", "end", "strand")

gene_df$strand[gene_df$strand == "+"] <- TRUE
gene_df$strand[gene_df$strand == "-"] <- FALSE
gene_df$`orf#` <- as.numeric(str_remove(gene_df$`orf#`, ">"))

gene_df$contig <- "KX687704.1"


# add annotation ----------------------------------------------------------

gene_df$annot <- "unknown"

gene_df$annot[7] <- "ORF6"
gene_df$annot[9] <- "ORF8"
gene_df$annot[11] <- "ORF10"
gene_df$annot[18] <- "ORF17"
gene_df$annot[22] <- "ORF21"
gene_df$annot[24] <- "ORF23"
gene_df$annot[25] <- "ORF24"
gene_df$annot[26] <- "ORF25"
gene_df$annot[8] <- "ORF7"


# plot genome map ---------------------------------------------------------

ggplot(gene_df, aes(xmin = start, xmax = end, y = contig, fill = annot)) +
    geom_gene_arrow() +
    facet_wrap(~ contig, scales = "free", ncol = 1) +
    theme_genes() +
    scale_fill_manual(values = c("ORF6" = "#f94144",
                                 "ORF8" = "#f3722c",
                                 "ORF10" = "#f8961e",
                                 "ORF17" = "#f9c74f",
                                 "ORF21" = "#90be6d",
                                 "ORF23" = "#43aa8b",
                                 "ORF24" = "#577590",
                                 "unknown" = "#FBFEF9",
                                 "ORF7" = "#D0D3D6",
                                 "ORF25" = "#D0D3D6")) +
    theme(legend.position = "None") +
    ylab("") +
    xlab("Position [bp]")

ggsave(plot = last_plot(), height = 4, width = 16, units = "cm",
       filename = "../plots/FIGURE_1_genome_map.svg")
