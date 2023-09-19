# author: dlueckin
# date: Tue Sep  5 12:14:57 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(googlesheets4)
library(gggenes)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# explore -----------------------------------------------------------------

df <- fread("FINAL_complete_annotations.csv") %>% 
    filter(Type == "repeat_region")


# clean up ----------------------------------------------------------------

# create unique identifier for each repat pair
df$repeat_id <- paste0(df$`Sequence Name`, "_", df$Name)

# number of reapats with the same name
df$repeat_count <- ave(df$repeat_id, df$repeat_id, FUN = length)

# remove singletons
df <- df %>% 
    filter(repeat_count >= 2)

# calculate distance for pairs
df <- df %>%
    group_by(repeat_id) %>%
    mutate(distance = ifelse(repeat_count == 2, abs(diff(Minimum)), NA))

# type of repeat (pair, multiple)
df <- df %>%
    mutate(repeat_type = case_when(
        repeat_count == 2 & distance == 0 ~ "palindromic",
        repeat_count == 2 & distance <= 200 ~ "pair_close",
        repeat_count == 2 & distance > 200 ~ "pair_distant",
        repeat_count > 2 ~ "multiple",
        TRUE ~ NA_character_))


# load region df ----------------------------------------------------------

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))
region_df <- region_df %>% 
    filter(contig_category == "complete") %>% 
    select(contig, sequence_length, region_start_bp, region_end_bp, region_direction)




# join dfs ----------------------------------------------------------------

names(df) <- c("name", "contig", "type", "start", "end", "length", "direction", "repeat_id", "repeat_count", "distance", "repeat_type")

big_df <- left_join(df, region_df, by = "contig")


big_df <- big_df %>%
    mutate(repeat_position = case_when(
        abs(start - region_start_bp) < abs(start - region_end_bp) & abs(start - region_start_bp) <= 20000 ~ "early",
        abs(start - region_end_bp) < abs(start - region_start_bp) & abs(start - region_end_bp) <= 20000 ~ "late",
        TRUE ~ "distant"
    ))

# flip early <-> late when region is reversed
big_df <- big_df %>%
    mutate(repeat_position = case_when(
        region_direction == "reverse" & repeat_position == "early" ~ "late",
        region_direction == "reverse" & repeat_position == "late" ~ "early",
        TRUE ~ repeat_position
    ))


# first we plot only the close repeats ------------------------------------
# set up mock db with a start and end (to trick gggenes into drawing the full thing)

genes_df <- region_df %>% 
    filter(!is.na(region_start_bp))

genes_df$annot <- "pR1SE_region"

tmp_df <- data.table()
for(i in 1:nrow(genes_df)){
    mock_start <- data.table("contig" = genes_df$contig[i],
                             "sequence_length" = 0,
                             "region_start_bp" = 0,
                             "region_end_bp" = 1,
                             "region_direction" = "forward",
                             "annot" = "start")
    mock_end <- data.table("contig" = genes_df$contig[i],
                           "sequence_length" = 0,
                           "region_start_bp" = genes_df$sequence_length[i] - 1,
                           "region_end_bp" = genes_df$sequence_length[i],
                           "region_direction" = "forward",
                           "annot" = "end")
    tmp_df <- rbind(tmp_df, rbind(mock_start, mock_end))
    
}

genes_df <- rbind(genes_df, tmp_df)
rm(tmp_df, mock_end, mock_start, df, i)


# clean the df ------------------------------------------------------------

genes_df <- genes_df %>% 
    select(-sequence_length)

genes_df$orientation <- ifelse(genes_df$direction == "forward", 1, 0)

names(genes_df) <- c("contig", "start", "end", "direction", "annotation",  "orientation")
    
    
    
# prepare repeats df for plotting -----------------------------------------

close_repeats_df <- big_df %>% 
    filter(repeat_type == "pair_close") %>% 
    group_by(repeat_id) %>% 
    sample_n(1) %>% 
    ungroup()

close_repeats_df <- close_repeats_df %>% 
    select(contig, repeat_position, start, end, direction)

close_repeats_df$orientation <- ifelse(close_repeats_df$direction == "forward", 1, 0)

names(close_repeats_df) <- c("contig", "annotation", "start", "end", "direction", "orientation")



# bind and plot -----------------------------------------------------------

plot_df <- rbind(genes_df, close_repeats_df)
plot_df$orientation <- ifelse(plot_df$direction == "forward", 1, 0)

# shorter names
plot_df$contig <- str_remove(plot_df$contig, pattern = "\\_length.*$")

# create dummies
dummies <- make_alignment_dummies(
    plot_df,
    aes(xmin = start, xmax = end, y = contig, id = annotation),
    on = "pR1SE_region"
)


# for 10 per plot
unique_contigs <- unique(plot_df$contig)

for(i in 1:4){
    current_contigs <- unique_contigs[(i*10-9):(i*10)]
    ggplot(plot_df %>% filter(contig %in% current_contigs), 
           aes(xmin = start, xmax = end+1000, y = contig, fill = annotation,
               forward = orientation)) +
        geom_gene_arrow(alpha = 0.8) +
        # geom_blank(data = dummies %>%  filter(contig %in% current_contigs)) +
        facet_wrap(~ contig, scales = "free", ncol = 1) +
        scale_fill_manual(values = c("early" = "blue",
                                     "late" = "orange",
                                     "start" = "white",
                                     "end" = "white",
                                     "distant" = "red",
                                     "pR1SE_region" = "darkgreen")) +
        theme_genes() +
        theme(legend.position = "None")
    
    ggsave(plot = last_plot(),
           file = paste0("../../plots/repeat_plots/plot_", i, ".pdf"), 
           height = 11, width = 7.7)
}



# plot distant vs contig_size ---------------------------------------------

repeat_df <- region_df
repeat_df <- left_join(repeat_df,
                       big_df %>% 
                           filter(repeat_type == "pair_distant") %>% 
                           group_by(contig) %>% 
                           count(),
                       by = "contig")
repeat_df$n <- repeat_df$n / 2
ggplot(repeat_df, aes(x = sequence_length, y = n)) +
    geom_point() +
    theme_cowplot() +
    xlim(c(0, 300000)) +
    ylim(c(0, 35)) +
    ylab("number of distant pairs")

ggsave("../../plots/number_of_distant_pairs_vs_sequence_length.png")                       


repeat_df <- left_join(repeat_df,
                       big_df %>% 
                           filter(repeat_type == "multiple") %>% 
                           group_by(contig) %>% 
                           count(),
                       by = "contig")
