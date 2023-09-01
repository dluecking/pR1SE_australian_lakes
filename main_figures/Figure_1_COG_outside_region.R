# author: dlueckin
# date: Wed Aug  9 19:14:32 2023

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


# load region_df ----------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(contig_category == "complete") %>% 
    filter(genome_segment == "circular") %>% 
    filter(integrated_state != "main_chromosome")



# load data ---------------------------------------------------------------

blast_df <- rbindlist(lapply(list.files("../J_COG_categories_complete/blast_out", full.names = TRUE), fread))


def <- fread("../../misc/test_new_COG/cog-20.def.tab")
# cog needs to be treated with sed -i 's/,,,,,/,,,,/' cog-20.cog.csv
cog <- fread("../../misc/test_new_COG/cog-20.cog.csv", fill = T, sep = ",") 

blast_df$COG <- cog$V7[match(blast_df$V2, cog$V3)]
blast_df$letter <- def$V2[match(blast_df$COG, def$V1)]

blast_df$orf <- as.numeric(str_extract(blast_df$V1, pattern = "\\d*$"))
blast_df$contig <- str_remove(blast_df$V1, pattern = "\\_\\d*$")


# keep only hits OUTSIDE of region borders --------------------------------

filtered_blast_df <- data.table()
for(i in 1:nrow(manual_region_df)){
    tmp_df <- blast_df %>% 
        filter(contig == manual_region_df$contig[i]) %>% 
        filter(!between(x = orf,
                        left = manual_region_df$manual_min_orf[i],
                        right = manual_region_df$manual_max_orf[i]))
    
    
    filtered_blast_df <- rbind(filtered_blast_df, tmp_df)
}
rm(tmp_df, i)

# add the counts in a template to plot ------------------------------------

to_plot <- fread("../../../template_classifier_count.tsv")

for(i in 1:nrow(to_plot)){
    current_letter <- to_plot$LETTER[i]
    to_plot$COUNT[i] <- filtered_blast_df %>% 
        filter(str_detect(string = letter, pattern = current_letter)) %>% 
        nrow()
}

# percentage
to_plot$perc <- to_plot$COUNT / sum(to_plot$COUNT)



# plot --------------------------------------------------------------------

ggplot(to_plot, aes(x = SHORT_DESCRIPTION, y = perc)) +
    geom_bar(stat = 'identity', position = "dodge", color = "black", fill = "#7F7F7F", alpha = 0.8) +
    xlab("COG Category") +
    ylab("Frequency") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45 ,  hjust=1),
          plot.margin=unit(c(0.5,0.5,0.5,l = 2), "cm"))

ggsave(plot = last_plot(), filename= "../plots/FIGURE_1_COG_barplot.svg", height = 4, width = 6)    
ggsave(plot = last_plot(), filename= "../plots/FIGURE_1_COG_barplot.png", height = 4, width = 6)  