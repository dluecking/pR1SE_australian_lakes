# author: dlueckin
# date: Wed Oct 11 13:27:52 2023

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



# import ------------------------------------------------------------------

df <- fread("plasmid_df_with_tax.csv")


# extract class -----------------------------------------------------------

df$class <- str_remove(str_remove(str_extract(df$taxonomy, "class\\_.*"), ",.*$"), "class\\_") 

df <- df%>% 
    filter(!is.na(class))


# new ---------------------------------------------------------------------
# I've got the number of genomes per class from GTDB, btw first one could also be 359 since there is thermoprotei and thermoprotei_A
division_factors <- data.frame(class = unique(df$class), division_factor = c(122, 750, 311, 
                                                                             141, 61, 559,
                                                                             801, 71)) 

plot_df <- as.data.table(table(df$class))
names(plot_df) <- c("class", "count")
plot_df <- left_join(plot_df, division_factors)

plot_df$plasmid_per_genome <- plot_df$count/plot_df$division_factor    


ggplot(plot_df, aes(x = class, y = plasmid_per_genome)) +
    geom_bar(stat = 'identity', color = "black", "fill" = "darkgrey", na.rm = TRUE) +
    coord_flip() +
    xlab("Class") +
    ylab("Plasmids / Genome") +
    theme_cowplot(
        font_size = 9
    )
    # ggtitle("Number of Plasmids per Class of Archaea")

ggsave("../../plots/number_of_plasmids_per_class_adjusted.png", plot = last_plot(), height = 2, width = 3)


# this is old, not adjusted per genome ------------------------------------

p1 <- ggplot(df, aes(x = class)) +
    geom_bar(color = "black", "fill" = "darkgrey", na.rm = TRUE) +
    coord_flip() +
    xlab("Class") +
    ylab("Count") +
    theme_cowplot() +
    ggtitle("Number of Plasmids per Class of Archaea")
# ggsave("../../plots/number_of_plasmids_per_class_total.png", plot = p1, height = 3, width = 8)

CUTOFF <- 200000
p2 <- ggplot(df %>% filter(length >= CUTOFF), aes(x = class)) +
    geom_bar(color = "black", "fill" = "darkgrey", na.rm = TRUE) +
    coord_flip() +
    xlab("Class") +
    ylab("Count") +
    theme_cowplot() +
    ggtitle("Number of Plasmids per Class of Archaea",
            subtitle = paste0("Length cutoff >=", CUTOFF))
# ggsave("../../plots/number_of_plasmids_per_class_cutoff.png", plot = p2, height = 3, width = 8)




