# author: dlueckin
# date: Tue Jun 27 11:59:11 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(googlesheets4)
library(cowplot)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# broder_df ---------------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(contig_category == "complete")


# read protein_df with annotations ----------------------------------------

prot_df <- fread("../G_annotate_complete_relatives/scripts/improved_gene_df.tsv")
prot_df <- prot_df %>% 
    filter(contig %in% manual_region_df$contig)

manual_region_df$insert_size <- 0
manual_region_df$inserted_genes <- 0

for(i in 1:nrow(manual_region_df)){
    ORF10_end <- prot_df %>% 
        filter(contig == manual_region_df$contig[i], annot == "ORF10.faa") %>% 
        select(end) %>% 
        pull()
    
    ORF17_start <- prot_df %>% 
        filter(contig == manual_region_df$contig[i], annot == "ORF17.faa") %>% 
        select(start) %>% 
        pull()

    manual_region_df$insert_size[i] <- abs(ORF17_start - ORF10_end)
    
    
    
    ORF10_orfn <- prot_df %>% 
        filter(contig == manual_region_df$contig[i], annot == "ORF10.faa") %>% 
        select(`orf#`) %>% 
        pull() %>% 
        as.numeric()
    
    ORF17_orfn <- prot_df %>% 
        filter(contig == manual_region_df$contig[i], annot == "ORF17.faa") %>% 
        select(`orf#`) %>% 
        pull() %>% 
        as.numeric()
    
    manual_region_df$inserted_genes[i] <- abs(ORF17_orfn - ORF10_orfn)
}


p1 <- ggplot(manual_region_df, aes(x = inserted_genes)) +
    geom_density() +
    xlim(c(0, 12)) +
    theme_cowplot(
        font_size = 10
    ) +
    xlab("Genes") +
    ylab("Density") +
    ggtitle("Number of genes inserted", 
            subtitle = paste0("complete apHPVs only (n = ", nrow(manual_region_df), ")"))


p2 <- ggplot(manual_region_df, aes(x = insert_size)) +
    geom_density() +
    xlim(c(0, 7500)) +
    theme_cowplot(
        font_size = 10
    ) +
    xlab("Length (bp)") +
    ylab("Density") +
    ggtitle("Insert size (bp)", 
            subtitle = paste0("complete apHPVs only (n = ", nrow(manual_region_df), ")"))
p3 <- p1 / p2

ggsave(plot = p3, filename = "../plots/insert_size.pdf", height = 6, width = 6)
ggsave(plot = p3, filename = "../plots/insert_size.png", height = 6, width = 6)
