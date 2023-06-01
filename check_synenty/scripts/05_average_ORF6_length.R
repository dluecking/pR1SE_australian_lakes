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

manual_region_df <- fread("manual_region_df.tsv") %>% 
    filter(comment == "")

gene_df <- fread("gene_df.tsv")
gene_df2 <- data.table()


for(i in 1:nrow(manual_region_df)){
    current_contig <- manual_region_df$contig[i]
    current_min <- manual_region_df$min_orf[i]
    current_max <- manual_region_df$max_orf[i]
    
    tmp_df <- gene_df %>% 
        filter(contig == current_contig) %>% 
        filter(`orf#` >= current_min, `orf#` <= current_max)
    
    # reverse order based on direction of ORF17
    if(tmp_df$strand[tmp_df$annot == "ORF17.faa"] != TRUE){
        tmp_df$end <- tmp_df$end * -1
        tmp_df$start <- tmp_df$start * -1
        
        # switch start and end to adjust arrow direction later on
        t <- tmp_df$start
        tmp_df$start <- tmp_df$end
        tmp_df$end <- t
    }
    gene_df2 <- rbind(gene_df2, tmp_df)
}
gene_df <- gene_df2
rm(gene_df2, i, current_contig, current_max, current_min, tmp_df)

# change ORF10 annot
gene_df$annot <- str_remove(gene_df$annot, "\\.faa")
gene_df$annot[gene_df$annot == "ORF10a"] <- "ORF10"
gene_df$annot[gene_df$annot == "ORF10b"] <- "ORF10"
gene_df$annot[gene_df$annot == ""] <- "unknown"



# calculate sizes of ORFs -------------------------------------------------

gene_df$length <- gene_df$end - gene_df$start


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
    xlab("Length (bp)") +
    ylab("")
ggsave(plot = box, "../../plots/average_ORF_lengths_ridgeline.pdf", height = 4, width = 7)



ridge <- ggplot(gene_df, aes(y = annot, x = length, fill = annot)) +
    geom_density_ridges(alpha=0.6, stat="binline", bins=70) +
    theme_ridges() +
    scale_fill_manual(values = c("ORF6" = "#f94144",
                                 "ORF8" = "#f3722c",
                                 "ORF10" = "#f8961e",
                                 "ORF17" = "#f9c74f",
                                 "ORF21" = "#90be6d",
                                 "ORF23" = "#43aa8b",
                                 "ORF24" = "#577590",
                                 "unknown" = "#FBFEF9")) +
    xlab("Length (bp)") +
    ylab("")


ggsave(plot = ridge, "../../plots/average_ORF_lengths_ridgeline.pdf", height = 4, width = 7)
