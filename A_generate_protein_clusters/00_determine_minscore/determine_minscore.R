# author: dlueckin
# date: Wed Jun 14 11:16:40 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggridges)

# TWO PARTS! 
# TOP: hmmsearch against proteins
# BOTTOM: hmmsearch against nr

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# load input --------------------------------------------------------------

hmm_df <- data.table()

for(f in list.files("hmm_out/", pattern = "hmm.domtblout$")){
    tmp_df <- fread(paste0("hmm_out/", f),skip = 3, fill = TRUE) %>% 
        select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
        filter(V1 != "#")
    names(tmp_df) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
    tmp_df$orf <- str_extract(string = f, pattern = "ORF\\d*")
    
    # calc alignment length
    tmp_df$aln_length <- tmp_df$alignment_end - tmp_df$alignment_start
    
    # prepare data for later filtering
    tmp_df$query_length <- as.numeric(tmp_df$query_length)
    tmp_df$score <- as.numeric(tmp_df$score)
    
    hmm_df <- rbind(hmm_df, tmp_df)
}
rm(f, tmp_df)



# plot score --------------------------------------------------------------

ggplot(hmm_df, aes(y = orf, x = score, fill = orf)) +
    geom_density_ridges() +
    theme_ridges() + 
    theme(
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
    ) +
    xlab("bitscore") +
    ylab("") +
    ggtitle("bitscore distribution", subtitle = "hmmsearch against uncurated blastp proteins")

ggsave(file = "../../plots/against_proteins_minscore_density.pdf", plot = last_plot(), height = 8, width = 16)

ggplot(hmm_df, aes(y = orf, x = score, fill = orf)) +
    geom_boxplot() +
    coord_flip() +
    theme_minimal()


# same but against nr ----------------------------------------------------

hmm_df <- data.table()

for(f in list.files("hmm_out/", pattern = "hmm.nr.domtblout")){
    
    is_not_empty <- as.numeric(str_split(system(command = paste0("wc -l hmm_out/", f), intern = TRUE), " ")[[1]][1]) > 13
    if(is_not_empty){
        tmp_df <- fread(paste0("hmm_out/", f),skip = 3, fill = TRUE) %>% 
            select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
            filter(V1 != "#")
        names(tmp_df) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
        tmp_df$orf <- str_extract(string = f, pattern = "ORF\\d*")
        
        # calc alignment length
        tmp_df$aln_length <- tmp_df$alignment_end - tmp_df$alignment_start
        
        # prepare data for later filtering
        tmp_df$query_length <- as.numeric(tmp_df$query_length)
        tmp_df$score <- as.numeric(tmp_df$score)
        
        hmm_df <- rbind(hmm_df, tmp_df)
    }
}
rm(f, tmp_df)

upper_bounds <- hmm_df %>% 
    group_by(orf) %>% 
    summarise(upper_bound = quantile(s, 0.75) + 1.5 * (quantile(score, 0.75) - quantile(score, 0.25) + 10))




ggplot(hmm_df, aes(y = orf, x = score, fill = orf)) +
    geom_density_ridges() +
    theme_ridges() + 
    # geom_point(aes(y = as.numeric(factor(c("ORF17", "ORF10"))) + 0.1, x = c(205, 120)), shape = 25, fill = "black") +
    theme(legend.position = "None")
    

