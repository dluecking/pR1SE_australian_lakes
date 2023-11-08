# author: dlueckin
# date: Tue Sep 19 17:50:15 2023

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


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
CUTOFF <- 5000


df <- region_df %>% 
    filter(integrated_state == "small_plasmid", contig_category == "complete")

# for me
blast_out <- rbindlist(lapply(list.files("blast_out", full.names = TRUE), fread))
names(blast_out) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", " "))

min95_2k <- blast_out %>% 
    filter(length >= CUTOFF) %>% 
    mutate(sseqid = str_remove(sseqid, "^NZ\\_")) %>% 
    filter(sseqid != qseqid)  %>% 
    distinct() %>% 
    filter(pident >= 80) %>% 
    filter(qseqid %in% df$contig)


# function ----------------------------------------------------------------

calculate_covered_length <- function(data, sequence_length, region_start, region_end) {
    v <- rep(0, sequence_length)
    for(i in 1:nrow(data)){
        v[data$qstart[i]:data$qend[i]] <- 1
    }
    # turn all elements within region to 0
    v[region_start:region_end] <- 0
    # count non-zero elements
    length_covered <- sum(v != 0)
}

# calculate ---------------------------------------------------------------

to_plot <- data.table()


for(contig in df$contig){
    region_s <- df$region_start_bp[df$contig == contig]
    region_e <- df$region_end_bp[df$contig == contig]
    sequence_length <- df$sequence_length[df$contig == contig]
    
    min95_2k <- blast_out %>% 
    filter(length >= CUTOFF) %>% 
        mutate(sseqid = str_remove(sseqid, "^NZ\\_")) %>% 
        filter(sseqid != qseqid)  %>% 
        distinct() %>% 
        filter(pident >= 80) %>% 
        filter(qseqid == contig) %>% 
        filter(!between(x = qstart, region_s, region_e)) %>% 
        filter(!between(x = qend, region_s, region_e))
    
    if(nrow(min95_2k) < 1){
        number_of_hits <- 0
        total_length <- 0
        first_five_sseqids <- ""
        
    }else{
        number_of_hits <- nrow(min95_2k)
        total_length <- calculate_covered_length(min95_2k, 
                                                 sequence_length = sequence_length, 
                                                 region_start = region_s, 
                                                 region_end = region_e)
        first_five_sseqids <- as.data.table(table(min95_2k$sseqid)) %>% 
            slice_max(order_by = N, n = 5) %>% 
            select(V1) %>% 
            pull()
    }
    print(paste0("Contig: ", contig, " - matching regions: ", number_of_hits, " - total length: ", total_length))
    print(first_five_sseqids, collapse = " -- ")
    
    
    tmp_df <- data.table(contig = contig,
                         region_length = df$region_length_bp[df$contig == contig],
                         sequence_length = sequence_length,
                         covered_length = total_length)
    to_plot <- rbind(to_plot, tmp_df)
    
}

to_plot$uncovered_length <- to_plot$sequence_length - to_plot$region_length - to_plot$covered_length

a <- gather(to_plot, key = "which", "length", -contig)
a$contig <- str_remove(a$contig, "\\_length.*")
# Ensure that 'which' is a factor variable
a$which <- as.factor(a$which)

# Create a new column for ordering by region_length_bp
a <- a %>%
    group_by(contig) %>%
    mutate(order = ifelse(which == "sequence_length", -length, length)) %>%
    ungroup()
# Order the contig variable
a$contig <- factor(a$contig, levels = unique(a$contig[order(a$order)]))

a <- a %>% filter(which != "sequence_length")

# plot --------------------------------------------------------------------


ggplot(a, aes(x = contig, y = length, fill = which)) +
    geom_bar(stat = 'identity', color = "black", alpha = 0.8) +
    coord_flip() +
    xlab("Contig") +
    ylab("Sequence Length (bp)") +
    labs(fill = "") +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("region_length" = "orange",
                                 "covered_length" = "#8338EC",
                                 "uncovered_length" = "grey"))









