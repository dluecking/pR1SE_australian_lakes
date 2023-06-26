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
library(googlesheets4)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# create automated region borders df --------------------------------------

gene_df <- fread("gene_df.tsv")

auto_df <- gene_df %>% 
    group_by(contig) %>% 
    summarise(min_orf = min(`orf#`), max_orf = max(`orf#`))
auto_df$region_length <- auto_df$max_orf - auto_df$min_orf
auto_df$manual_min_orf <- NA
auto_df$manual_max_orf <- NA
auto_df$manual_region_length <- NA

fwrite(auto_df, "automated_region_borders.tsv")

# I ASSUME YOU HAVE SETUP THE MANUAL_REGION_DF BASED ON THE AUTOMATED ONE BUT WITH MANUAL CURATION
rm(auto_df)



# import df and filter  --------------------------------------------------


manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                               sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(visualize == "TRUE")


# categorize contigs ------------------------------------------------------

manual_region_df$number_of_unique_hits  <- 0
for(i in 1:nrow(manual_region_df)){
    manual_region_df$number_of_unique_hits[i] <-  gene_df %>% 
        filter(contig == manual_region_df$contig[i]) %>% 
        summarise(length(unique(annot))) %>% 
        pull() - 1
}

manual_region_df$contig_category <- "fragment"
manual_region_df$contig_category[manual_region_df$number_of_unique_hits >= 6] <- "complete"


# plot fragment ones ------------------------------------------------------

gene_df_complete <- data.table()
manual_complete <- manual_region_df[manual_region_df$contig_category == "complete"]


for(i in 1:nrow(manual_complete)){
    current_contig <- manual_complete$contig[i]
    current_min <- manual_complete$manual_min_orf[i]
    current_max <- manual_complete$manual_max_orf[i]
    
    tmp_df <- gene_df %>% 
        filter(contig == current_contig) %>% 
        filter(`orf#` >= current_min, `orf#` <= current_max)
    
    # reverse order based on direction of ORF21
    if(tmp_df$strand[tmp_df$annot == "ORF21.faa"] != TRUE){
        tmp_df$end <- tmp_df$end * -1
        tmp_df$start <- tmp_df$start * -1

        # switch start and end to adjust arrow direction later on
        t <- tmp_df$start
        tmp_df$start <- tmp_df$end
        tmp_df$end <- t
    }
    gene_df_complete <- rbind(gene_df_complete, tmp_df)
}

rm(i, current_contig, current_max, current_min, tmp_df)


# plot --------------------------------------------------------------------

gene_df_complete$annot <- str_remove(gene_df_complete$annot, "\\.faa")
gene_df_complete$annot[gene_df_complete$annot == ""] <- "unknown"

# remove long ass names
gene_df_complete$contig <- str_remove(gene_df_complete$contig, "\\_length.*")
gene_df_complete$length <- gene_df_complete$end - gene_df_complete$start

# save that for later
fwrite(gene_df_complete, "cleaned_gene_df_complete.tsv")


# create dummies
dummies <- make_alignment_dummies(
    gene_df_complete,
    aes(xmin = start, xmax = end, y = contig, id = annot),
    on = "ORF24"
)

# plotting in 4 parts
for(i in 1:5){
    current_contigs <- unique(gene_df_complete$contig)[(i*8-7):(i*8)]
    current_plot <- ggplot(gene_df_complete %>% 
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
    
    ggsave(current_plot, filename = paste0("../../plots/gene_map_plot_complete", i, ".pdf"), height = 7.7, width = 11)
}
rm(dummies, current_plot, current_contigs, i, t)


# plot fragments ----------------------------------------------------------

gene_df_fragment <- data.table()
manual_fragment <- manual_region_df[manual_region_df$contig_category == "fragment"]

for(i in 1:nrow(manual_fragment)){
    current_contig <- manual_fragment$contig[i]
    current_min <- manual_fragment$manual_min_orf[i]
    current_max <- manual_fragment$manual_max_orf[i]
    
    tmp_df <- gene_df %>% 
        filter(contig == current_contig) %>% 
        filter(`orf#` >= current_min, `orf#` <= current_max)
    
    gene_df_fragment <- rbind(gene_df_fragment, tmp_df)
}

rm(i, current_contig, current_max, current_min, tmp_df)


# cosmetics ---------------------------------------------------------------

gene_df_fragment$annot <- str_remove(gene_df_fragment$annot, "\\.faa")
gene_df_fragment$annot[gene_df_fragment$annot == ""] <- "unknown"

# remove long ass names
gene_df_fragment$contig <- str_remove(gene_df_fragment$contig, "\\_length.*")
gene_df_fragment$length <- gene_df_fragment$end - gene_df_fragment$start

# save that for later
fwrite(gene_df_fragment, "cleaned_gene_df_fragment.tsv")


# plot --------------------------------------------------------------------

# plotting in 8 parts
for(i in 1:4){
    current_contigs <- unique(gene_df_fragment$contig)[(i*8-7):min((i*8), nrow(manual_fragment))]
    current_plot <- ggplot(gene_df_fragment %>% 
                               filter(contig %in% current_contigs), aes(xmin = start, xmax = end, y = contig, fill = annot)) +
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
                                     "unknown" = "#FBFEF9"))
    
    ggsave(current_plot, filename = paste0("../../plots/gene_map_plot_fragment_", i, ".pdf"), height = 7.7, width = 11)
}
rm(current_plot, current_contigs, i)
