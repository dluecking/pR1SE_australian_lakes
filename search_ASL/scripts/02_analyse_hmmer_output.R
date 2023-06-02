# author: dlueckin
# date: Tue May  9 13:14:32 2023

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

core_ORFs <- c("ORF6", "ORF24", "ORF23", "ORF21")


# explore -----------------------------------------------------------------

lakes <- fread("list_of_lakes.txt", header = FALSE)
orf_df <- fread("ORF_cutoff_table.tsv")
big_df <- data.table()

for(file in list.files("../01_hmm_results", full.names = TRUE, pattern = "domtblout.txt")){
    is_not_empty <- as.numeric(str_split(system(command = paste0("wc -l ", file), intern = TRUE), " ")[[1]][1]) > 13
    
    if(is_not_empty){
        domtbl <- fread(file, skip = 3, fill = TRUE) %>% 
            select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
            filter(V1 != "#")
        
        names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
        
        current_lake <- str_extract(file, pattern = paste(lakes$V1, collapse = "|"))
        domtbl$lake <- current_lake
        current_ORF <- str_extract(file, pattern = "ORF[^\\_]*")
        
        # calc alignment length
        domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start
        
        # prepare data for later filtering
        domtbl$query_length <- as.numeric(domtbl$query_length)
        domtbl$score <- as.numeric(domtbl$score)
        
        domtbl <- domtbl %>% 
            filter(score >= 50)
            # filter(score >= orf_df[orf_df$orf == current_ORF]$min_score) %>% 
            # filter(aln_length >= orf_df[orf_df$orf == current_ORF]$min_alignment_length)
        
        big_df <- rbind(big_df, domtbl)
    }
}

<<<<<<< HEAD:search_ASL/scripts/analyse_hmmer_output.R
rm(orf_df, lakes, domtbl, current_lake, current_ORF)



# create contig_df --------------------------------------------------------

contig_df <- data.table(table(big_df$hit))
contig_df$contig <- str_remove(contig_df$V1, pattern = "\\_length.*")
contig_df$contig_len <- as.numeric(str_remove(str_extract(contig_df$V1, "length\\_\\d*"), "length\\_"))

contig_df <- contig_df %>% 
    filter(N >= 4)
=======
big_df$contig <- str_remove(big_df$hit, pattern = "\\_length.*")
big_df$ORF <- str_remove(big_df$query, "\\.faa")

rm(file, is_not_empty, current_lake, current_ORF, orf_df, lakes, domtbl)


# generate presence absence plot ------------------------------------------

wide_df <- data.table(table(big_df$contig))
wide_df <- wide_df %>% 
    filter(N >= 2)
names(wide_df) <- c("contig", "N")

big_df$presence <- 1
# Cast the data to a wide format
wide_df <- dcast(big_df %>% filter(contig %in% wide_df$contig), contig ~ ORF, value.var = "presence", fill = 0)

wide_df[,2:8] <- lapply(wide_df[,2:8], function(x) ifelse(x > 1, 1, x))
wide_df$hits <- rowSums(wide_df[,2:8])
wide_df <- wide_df %>% 
    arrange(desc(hits)) %>% 
    filter(hits >= 3)

ggplot(melt(wide_df) %>% filter(variable != "hits"), aes(variable, contig, fill = variable, alpha = value)) + 
    geom_tile(colour = "gray50") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() + 
    coord_flip() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("presence_absence_plot_novel.png", plot = last_plot(), height = 8, width = 16)
fwrite(wide_df, "novel_hosts_df.tsv")





