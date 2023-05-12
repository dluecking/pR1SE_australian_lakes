# author: dlueckin
# date: Wed May 10 14:11:11 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(reshape2)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# import df and spread ----------------------------------------------------

prot_df <- fread("prot_df.tsv")

wide_df <- data.table(organism = unique(prot_df$organism))

prot_df$presence <- 1


# Cast the data to a wide format
wide_df <- dcast(prot_df, organism ~ orf, value.var = "presence", fill = 0)
wide_df[,2:9] <- lapply(wide_df[,2:9], function(x) ifelse(x > 1, 1, x))

wide_df$hits <- rowSums(wide_df[,2:9])
wide_df <- wide_df %>% 
    arrange(desc(hits)) %>% 
    filter(hits >= 4)


ggplot(melt(wide_df) %>% filter(variable != "hits"), aes(variable, organism, fill = variable, alpha = value)) + 
    geom_tile(colour = "gray50") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() + 
    coord_flip() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("presence_absence_plot.png", plot = last_plot(), height = 8, width = 16)


fwrite(wide_df, "known_hosts_df.tsv")
