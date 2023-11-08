# author: dlueckin
# date: Tue Oct 24 14:28:44 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# -------------------------------------------------------------------------

dirs <- list.dirs(".")
dirs <- dirs[str_detect(dirs, "ORF")]
dirs <- dirs[!str_detect(dirs, "protein")]
dirs <- dirs[!str_detect(dirs, "25")]

df <- data.table("files" = as.character(),
                 "orf" = as.character(),
                 "pdb_score" = as.numeric())

for(dir in dirs){
    files <- list.files(dir, pattern = "pdb")
    
    tmp <- data.table("files" = files,
                      "orf" = dir,
                      "pdb_score" = 0)
    
    df <- rbind(df, tmp)

}

df$full_path <- paste0(df$orf, "/", files)



# calculate score ---------------------------------------------------------

for(i in 1:nrow(df)){
    pdb <- fread(df$full_path[i], skip = 22)
    df$pdb_score[i] <- mean(pdb$V11)

}



# plot pdb scores ---------------------------------------------------------

ggplot(df, aes(x = orf, y = pdb_score)) + 
    geom_boxplot() +
    theme_cowplot() +
    ylab("mean confidence score [0-1]") +
    xlab("ORF")

ggsave(plot = last_plot(), file = "../../plots/pdb_scores_of_ORF_structures.png", height = 4, width = 5)
