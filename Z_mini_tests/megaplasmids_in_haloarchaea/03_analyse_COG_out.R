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


# load blast_df -----------------------------------------------------------

blast_df <- rbindlist(lapply(list.files("blast_out", full.names = TRUE), fread))


# add COG db, based on hit ------------------------------------------------

def <- fread("../../../misc/mvome_cog_analysis/helper_files/cog-20.def.tab")
# cog needs to be treated with sed -i 's/,,,,,/,,,,/' cog-20.cog.csv
cog <- fread("../../../misc/mvome_cog_analysis/helper_files/cog-20.cog.csv", fill = T, sep = ",") 

blast_df$COG <- cog$V7[match(blast_df$V2, cog$V3)]
blast_df$letter <- def$V2[match(blast_df$COG, def$V1)]


# filter blast_df by 10-5  ------------------------------------------------

blast_df <- blast_df %>% 
    filter(V11 <= 10^-5)


# prepare for plotting ----------------------------------------------------

template_df <- fread("../../../../template_classifier_count.tsv")

for(i in 1:nrow(template_df)){
    current_letter <- template_df$LETTER[i]
    template_df$COUNT[i] <- blast_df %>% 
        filter(str_detect(string = letter, pattern = current_letter)) %>% 
        nrow() %>% 
        unlist()
}
template_df$perc <- template_df$COUNT / sum(template_df$COUNT) * 100

fwrite(template_df, "df_to_plot_for_COG.tsv")
template_df <- fread("df_to_plot_for_COG.tsv")
# plot
ggplot(template_df, aes(x = reorder(SHORT_DESCRIPTION, desc(perc)), y = perc)) +
    geom_bar(stat = 'identity', color = "black", fill = "grey", alpha = 0.9) +
    theme_minimal() +
    xlab("COG category") +
    ylab("Assigned genes [%]") +
    ggtitle("COG category of genes on megaplasmids",
            subtitle = "5 plasmids") +
    theme(
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45 ,  hjust=1),
          plot.margin=unit(c(0.5,0.5,0.5,l = 3.5), "cm"))

ggsave(plot = last_plot(), filename = "../../plots/COG_megaplasmids.png", height = 4, width = 7,
       bg = "transparent")
ggsave(plot = last_plot(), filename = "../../plots/COG_megaplasmids.svg", height = 4, width = 7,
       bg = "transparent")
ggsave(plot = last_plot(), filename = "../../plots/COG_megaplasmids.pdf", height = 4, width = 7,
       bg = "transparent")

