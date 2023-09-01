# author: dlueckin
# date: Tue Aug 15 10:55:11 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(googlesheets4)
library(cowplot)


# working directory -------------------------------------------------------

this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# load region-df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))


# import interpro ---------------------------------------------------------

interpro_df <- fread("interpro_out_full.tsv")


# set up contig_df and filter ---------------------------------------------

contig_df <- data.table("contig" = unique(interpro_df$contig))

SEARCH_TERMS <- c("ParAB", "ParB", "ParC", "MinD", "Tubulin", "FtsZ")
# search for key words
for(word in SEARCH_TERMS){
    tmp_df <- interpro_df %>% 
        group_by(contig) %>% 
        summarize(word = ifelse(
            any(str_detect(annotation, regex(word, ignore_case = TRUE))
                ), TRUE, FALSE))
    names(tmp_df) <- c("contig", word)
    contig_df <- left_join(contig_df, tmp_df, by = "contig")
}

contig_df <- contig_df %>% 
    filter(contig %in% region_df$contig[region_df$contig_category == "complete"])


# combine information
to_plot <- gather(contig_df, key = "gene", value = "val", -contig)


# plot
p1 <- ggplot(to_plot, aes(x = gene,y = contig, alpha = val)) + 
    geom_tile(colour = "black", fill = "blue") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "None") +
    ggtitle("Partition Systems") +
    xlab("")

# ggsave("../../plots/partition_systems.pdf", height = 10, width = 10)




# same but for toxin antitoxin --------------------------------------------

SEARCH_TERMS <- c("RelE-like", "ParE", "AbrB", "MazE", "MazeF", "VapB", "PhoU", "PemI", "SpoVT", "toxin", "antitoxin")
contig_df <- data.table("contig" = unique(interpro_df$contig))

for(word in SEARCH_TERMS){
    tmp_df <- interpro_df %>% 
        group_by(contig) %>% 
        summarize(word = ifelse(
            any(str_detect(annotation, regex(paste0(".*", word, ".*"), ignore_case = TRUE))
            ), TRUE, FALSE))
    names(tmp_df) <- c("contig", word)
    contig_df <- left_join(contig_df, tmp_df, by = "contig")
}

contig_df <- contig_df %>% 
    filter(contig %in% region_df$contig[region_df$contig_category == "complete"])


# combine information
to_plot <- gather(contig_df, key = "gene", value = "val", -contig)

# plot
p2 <- ggplot(to_plot, aes(x = gene,y = contig, alpha = val)) + 
    geom_tile(colour = "black", fill = "red") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          legend.position = "None") +
    ggtitle("Toxin/Antitoxin Systems") +
    ylab("")

# ggsave("../../plots/toxin-antitoxin_systems.pdf", height = 10, width = 10)




# others ------------------------------------------------------------------

SEARCH_TERMS <- c("helicase", "cdc6", "orc1", "HNH", "endonuclease", "recombinase", "integrase", "transposase", "pilin", "Vaccinia Virus protein VP39")
contig_df <- data.table("contig" = unique(interpro_df$contig))

for(word in SEARCH_TERMS){
    tmp_df <- interpro_df %>% 
        group_by(contig) %>% 
        summarize(word = ifelse(
            any(str_detect(annotation, regex(paste0(".*", word, ".*"), ignore_case = TRUE))
            ), TRUE, FALSE))
    names(tmp_df) <- c("contig", word)
    contig_df <- left_join(contig_df, tmp_df, by = "contig")
}

contig_df <- contig_df %>% 
    filter(contig %in% region_df$contig[region_df$contig_category == "complete"])

# combine information
to_plot <- gather(contig_df, key = "gene", value = "val", -contig)

# plot
p3 <- ggplot(to_plot, aes(x = gene,y = contig, fill = gene, alpha = val)) + 
    geom_tile(colour = "black") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          legend.position = "None") +
    ggtitle("Others") +
    ylab("") +
    xlab("")
    

# ggsave("../../plots/other_genes_presence_absence.pdf", height = 10, width = 10)

# transposase stuff -------------------------------------------------------

SEARCH_TERMS <- c("Putative transposase DNA-binding domain", "Transposase DDE domain", "Probable transposase", 
                  "IS200", "ISH3", "IS66", "IS5", "IS4")
contig_df <- data.table("contig" = unique(interpro_df$contig))

for(word in SEARCH_TERMS){
    tmp_df <- interpro_df %>% 
        group_by(contig) %>% 
        summarize(word = ifelse(
            any(str_detect(annotation, regex(paste0(".*", word, ".*"), ignore_case = TRUE))
            ), TRUE, FALSE))
    names(tmp_df) <- c("contig", word)
    contig_df <- left_join(contig_df, tmp_df, by = "contig")
}

contig_df <- contig_df %>% 
    filter(contig %in% region_df$contig[region_df$contig_category == "complete"])


# combine information
to_plot <- gather(contig_df, key = "gene", value = "val", -contig)

# plot
p4 <- ggplot(to_plot, aes(x = gene,y = contig, alpha = val)) + 
    geom_tile(colour = "black", fill = "darkgreen") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          legend.position = "None") +
    ggtitle("Others") +
    ylab("") +
    xlab("")



p1 + p2 + p3

ggsave(plot = last_plot(), height = 10, width = 20, file = "../../plots/presence_absence_toxin_partition_others.pdf")
