# author: dlueckin
# date: Tue Jun 27 16:38:47 2023

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



# load important dfs ------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(contig_category == "complete")


# load blast_df -----------------------------------------------------------

blast_df <- rbindlist(lapply(list.files("blast_out", full.names = TRUE), fread))


# add COG db, based on hit ------------------------------------------------

def <- fread("../../misc/test_new_COG/cog-20.def.tab")
# cog needs to be treated with sed -i 's/,,,,,/,,,,/' cog-20.cog.csv
cog <- fread("../../misc/test_new_COG/cog-20.cog.csv", fill = T, sep = ",") 

blast_df$COG <- cog$V7[match(blast_df$V2, cog$V3)]


# add actual letter per hit -----------------------------------------------

blast_df$letter <- def$V2[match(blast_df$COG, def$V1)]


# add orf number, extract contig, keep only region genes ------------------

blast_df$orf <- as.numeric(str_extract(blast_df$V1, pattern = "\\d*$"))
blast_df$contig <- str_remove(blast_df$V1, pattern = "\\_\\d*$")

filtered_blast_df <- data.table()
for(i in 1:nrow(manual_region_df)){
    tmp_df <- blast_df %>% 
        filter(contig == manual_region_df$contig[i]) %>% 
        filter(orf >= manual_region_df$manual_min_orf[i],
               orf <= manual_region_df$manual_max_orf[i])
    
    filtered_blast_df <- rbind(filtered_blast_df, tmp_df)
}


to_plot <- as.data.table(table(filtered_blast_df$letter))

ggplot(to_plot, aes(x = V1, y = N)) +
    geom_bar(stat = 'identity') +
    theme_minimal() +
    ggtitle("COG category of genes on pR1SE relatives",
            subtitle = "202 knowns, 901 unknowns, 41 complete relatives")
ggsave(plot = last_plot(), filename = "../plots/COG_analysis.pdf", height = 5, width = 5)

