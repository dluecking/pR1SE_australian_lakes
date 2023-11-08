# author: dlueckin
# date: Tue Jun 27 16:38:47 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)
library(cowplot)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# load important dfs ------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))

VIABLE_CONTIGS <- region_df %>% 
    filter(contig_category == "complete") %>% 
    filter(!is.na(manual_min_orf)) %>% 
    filter(integrated_state != "main_chromosome") %>% 
    select(contig) %>% unlist()



# load COG dfs ------------------------------------------------------------

def <- fread("../../misc/mvome_cog_analysis/helper_files/cog-20.def.tab")
# cog needs to be treated with sed -i 's/,,,,,/,,,,/' cog-20.cog.csv
cog <- fread("../../misc/mvome_cog_analysis/helper_files/cog-20.cog.csv", fill = T, sep = ",") 


# load blast out ----------------------------------------------------------
# pR1SE out
blast_df <- rbindlist(lapply(list.files("blast_out_pR1SE/", full.names = TRUE), fread))
names(blast_df) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", " "))

blast_df$orf <- as.numeric(str_extract(blast_df$qseqid, pattern = "\\d*$"))
blast_df$contig <- str_remove(blast_df$qseqid, pattern = "\\_\\d*$")

# keep only viable contigs
blast_df <- blast_df %>% 
    filter(contig %in% VIABLE_CONTIGS)

# keep only outside of region
blast_df <- blast_df %>% 
    filter(!between(x = orf,
                    left = region_df$manual_min_orf[match(blast_df$contig, region_df$contig)] + 1,
                    right = region_df$manual_max_orf[match(blast_df$contig, region_df$contig)] - 1))
blast_df$origin <- "pR1SE"


# host out
tmp_df <- fread("blast_out_archaeal_proteins/blast_vs_COG.out")
names(tmp_df) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", " "))

tmp_df$orf <- 0
tmp_df$contig <- tmp_df$qseqid

tmp_df$origin <- "host"

# combine
blast_df <- rbind(blast_df, tmp_df)
rm(tmp_df)

# filter bad hits
blast_df <- blast_df %>% 
    filter(evalue <= 0.00001)


# combine COG info with blast df ------------------------------------------

blast_df$COG <- cog$V7[match(blast_df$sseqid, cog$V3)]
blast_df$letter <- def$V2[match(blast_df$COG, def$V1)]


# create COG plot ---------------------------------------------------------

# we do this, so the counting later works faster
count_df <- blast_df %>% 
    select(origin, letter) %>% 
    group_by(origin, letter) %>% 
    tally()


# load and repeat the template df
COG_template_df <- fread("../../../template_classifier_count.tsv")
repeat_df <- function(d, n) {
    return(do.call("rbind", replicate(n, d, simplify = FALSE)))
}
COG_template_df <- repeat_df(COG_template_df, length(unique(blast_df$origin))) 
COG_template_df$label <- rep(unique(blast_df$origin), each = nrow(COG_template_df) / length(unique(blast_df$origin)))

# go over each category and count
for(i in 1:nrow(COG_template_df)){
    COG_template_df$COUNT[i] <- count_df %>% 
        filter(str_detect(string = origin, pattern = COG_template_df$label[i])) %>%     # matches origin
        filter(str_detect(string = letter, pattern = COG_template_df$LETTER[i])) %>%    # matches label
        ungroup() %>%                                                                   # ungroup, otherwise the origin appears in the select
        select(n) %>% 
        summarise(total = sum(n)) %>% 
        pull()
}


# percentage wise
COG_template_df <- COG_template_df %>%
    group_by(label) %>%
    mutate(total_count = sum(COUNT)) %>%
    ungroup() %>%
    mutate(percentage = COUNT / total_count * 100)


# save df for plotting ----------------------------------------------------

fwrite(COG_template_df, "df_for_plotting.csv")


