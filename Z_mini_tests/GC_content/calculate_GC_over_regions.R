# author: dlueckin
# date: Wed Jul 26 14:35:54 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives")) %>% 
    filter(integrated_state == "small_plasmid", contig_category == "complete")



# plot gc over region -----------------------------------------------------

get_gc_in_window <- function(pos, sequence, WINDOW_SIZE){
    start <- max(pos-WINDOW_SIZE/2, 0)
    end <- min(pos + WINDOW_SIZE/2, length(sequence[[1]]))
    
    return(GC(sequence[[1]][start:end]))
    
}

region_df$gc_in_region <- 0

for(i in 1:nrow(region_df)){
    if(is.na(start <- region_df$region_start_bp[i]))
        next
    
    contig <- region_df$contig[i]
    start <- region_df$region_start_bp[i]
    end <- region_df$region_end_bp[i]
    
    sequence_path <- paste0("../../FINAL_complete_relatives/", contig, ".fasta")
    fasta <- read.fasta(sequence_path)
    sequence <- getSequence(fasta)
    
    region_df$gc_in_region[i] <- GC(sequence[[1]][start:end]) * 100
    
    gc_df <- data.table(pos = seq(1, length(sequence[[1]]), 100),
                        gc = 0)
    
    for(j in 1:nrow(gc_df)){
        gc_df$gc[j] <- get_gc_in_window(gc_df$pos[j], sequence, 100)
    }
    
    
    p1 <- ggplot(data = gc_df, aes(x = pos, y = gc)) +
        geom_line() +
        annotate("rect", xmin = start, xmax = end, ymin = 0, ymax = 1, alpha = 0.1, fill = "red") +
        ggtitle(paste0(contig, " (", region_df$genome_segment[i], ")"), 
                subtitle = paste0("GC = ", region_df$`% GC`[i], " - region_GC = ", region_df$gc_in_region[i])) +
        theme_minimal()
    ggsave(filename = paste0("../../plots/GC_plots/", contig, "_gc.pdf"), plot = p1, height = 3, width = 7)
    ggsave(filename = paste0("../../plots/GC_plots/", contig, "_gc.png"), plot = p1, height = 3, width = 7)
}
