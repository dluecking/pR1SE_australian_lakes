# author: dlueckin
# date: Mon Aug  7 14:36:39 2023

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



# load pilup data ---------------------------------------------------------

pile_df <- rbindlist(lapply(list.files("pileup/", full.names = TRUE), fread))
names(pile_df) <- c("seq", "pos", "ref", "depth", "bases", "qual")
pile_df$errors <- str_count(pile_df$bases, pattern = "A|C|T|G|a|c|t|g|\\*")
pile_df$correct_percentage <- (pile_df$depth - pile_df$errors) / pile_df$depth * 100
pile_df$correct_percentage[pile_df$correct_percentage < 0] <- 0



# load prot_df ------------------------------------------------------------

prot_df <-fread("../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv") %>% 
    filter(contig %in% str_remove(list.files("contigs/"), "\\.fasta")) %>% 
    select(contig, annot, start, end) %>% 
    filter(annot != "")


# calculate average correctness per gene ----------------------------------

prot_df$mean_correct <- 0
prot_df$positions_under95_correct <- 0

for(i in 1:nrow(prot_df)){
    tmp_df <- pile_df %>% 
        filter(seq == prot_df$contig[i]) %>% 
        filter(pos >= prot_df$start[i], pos <= prot_df$end[i]) %>% 
        select(correct_percentage)
    
    prot_df$mean_correct[i] <- mean(tmp_df$correct_percentage)
    prot_df$positions_under95_correct[i] <- tmp_df %>% filter(correct_percentage <= 95) %>% nrow() / nrow(tmp_df)
}


# plot data ---------------------------------------------------------------

p1 <- ggplot(prot_df, aes(x = annot, y = mean_correct)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(prot_df, aes(x = annot, y = positions_under95_correct)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1 + p2
ggsave("../plots/per_position_conservation_per_gene.png", height = 5, width = 5, plot = last_plot())




# calc bin df -------------------------------------------------------------

WINDOW <- 33
bin_df <- data.table()

for(i in 1:nrow(prot_df)){
    tmp_df <- pile_df %>% 
        filter(seq == prot_df$contig[i]) %>% 
        filter(pos >= prot_df$start[i], pos <= prot_df$end[i]) %>% 
        select(pos, correct_percentage)
    
    tmp_df$qual[tmp_df$correct_percentage <= 98] <- "<98"
    
    
    
    to_plot <- data.table(bin = 1:round(nrow(tmp_df)/WINDOW + 1, 0),
                          no_of_under_98 = 0)
    
    for(j in 1:nrow(to_plot)){
        current_bin <- to_plot$bin[j]
        
        to_plot$no_of_under_98[j] <- 
            tmp_df[(WINDOW*current_bin-WINDOW):(current_bin*WINDOW)] %>% 
            filter(!is.na(qual)) %>% 
            nrow()
        
        
    }
    
    to_plot$annot <- prot_df$annot[i]
    to_plot$contig <- prot_df$contig[i]
    
    bin_df <- rbind(bin_df, to_plot)
    
}


# plot bin df -------------------------------------------------------------

bin_df$contig <- str_remove(bin_df$contig, "\\_length.*$")
ggplot(bin_df , aes(x = bin, y = no_of_under_98, fill = contig)) +
    # geom_point(alpha = 0.7) +
    geom_bar(stat = 'identity', alpha = 0.7, color = "black") +
    facet_wrap(~annot, ncol = 1, strip.position = "right", scales = "free") + # scales free changes interpretation imo
    theme_minimal() +
    theme(legend.position = "bottom") +
    xlab(paste0("#Bin (bin size = ",WINDOW, " bp)")) +
    ggtitle(label = "Conserved regions along genes", subtitle = paste0("Number of positions within a ", WINDOW, " bp bin that have \nless than 98 % ID across mapped reads."))

ggsave(last_plot(), filename = "../plots/conserved_genes_bin_plot.png", height = 10, width = 9)
