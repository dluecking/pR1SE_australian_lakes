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


# plot --------------------------------------------------------------------

# change ORF10 annot
gene_df$annot <- str_remove(gene_df$annot, "\\.faa")
gene_df$annot[gene_df$annot == "ORF10a"] <- "ORF10"
gene_df$annot[gene_df$annot == "ORF10b"] <- "ORF10"
gene_df$annot[gene_df$annot == ""] <- "unknown"

# remove long ass names
gene_df$contig <- str_remove(gene_df$contig, "\\_length.*")


# create dummies
dummies <- make_alignment_dummies(
    gene_df,
    aes(xmin = start, xmax = end, y = contig, id = annot),
    on = "ORF24"
)

# plotting in 4 parts
for(i in 1:4){
    current_contigs <- unique(gene_df$contig)[(i*8-7):(i*8)]
    current_plot <- ggplot(gene_df %>% 
                               filter(contig %in% current_contigs), aes(xmin = start, xmax = end, y = contig, fill = annot)) +
        geom_gene_arrow() +
        geom_blank(data = dummies %>%  filter(contig %in% current_contigs)) +
        facet_wrap(~ contig, scales = "free", ncol = 1) +
        theme_genes() +
        scale_fill_manual(values = c("ORF6" = "#f94144",
                                     "ORF8" = "#f3722c",
                                     "ORF10" = "#f8961e",
                                     "ORF17" = "#f9c74f",
                                     "ORF21" = "#90be6d",
                                     "ORF23" = "#43aa8b",
                                     "ORF24" = "#577590",
                                     "unknown" = "#FBFEF9"))
    
    ggsave(current_plot, filename = paste0("plot", i, ".pdf"), height = 7.7, width = 11)
}






# old --------------------------------------------------------------------

h_df <- data.table("contig" = unique(gene_df$contig))
h_df$min_orf <- 0
h_df$max_orf <- 0

for(i in 1:nrow(h_df)){
    h_df$min_orf[i] <- min(gene_df$`orf#`[gene_df$contig == h_df$contig[i]])
    h_df$max_orf[i] <- max(gene_df$`orf#`[gene_df$contig == h_df$contig[i]])
}
fwrite(h_df, "automated_region_borders.tsv")




