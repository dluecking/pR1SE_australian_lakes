# author: dlueckin
# date: Thu Jun  1 10:31:21 2023

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



# load the two dfs --------------------------------------------------------

gair107 <- fread("torrens_mapped_to_Gairdner_107_basecov.txt")
gair107$contig <- "Gair107"
gair167 <- fread("torrens_mapped_to_Gairdner_167_basecov.txt")
gair167$contig <- "Gair167"
torrens <- fread("torrens_mapped_to_torrens_basecov.txt")
torrens$contig <- "torrens"

big_df <- rbind(gair107, gair167, torrens)



# plot --------------------------------------------------------------------

ggplot(big_df %>% filter(contig == "torrens"), aes(x = Pos, y = Coverage, group = contig, color = contig)) +
    geom_line()
