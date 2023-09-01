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

prot_df <- fread("prot_df.tsv") %>% 
    filter(organism_name != "")

wide_df <- data.table(organism = unique(prot_df$organism_acc))

prot_df$presence <- 1


# Cast the data to a wide format
wide_df <- dcast(prot_df, organism_acc ~ orf, value.var = "presence", fill = 0)
wide_df[,2:8] <- lapply(wide_df[,2:8], function(x) ifelse(x > 1, 1, x))

wide_df$hits <- rowSums(wide_df[,2:8])
wide_df <- wide_df %>% 
    arrange(desc(hits)) %>% 
    filter(hits >= 2)

# convert back to long
to_plot <- melt(wide_df) %>% filter(variable != "hits")

# add name 
to_plot$organism_name <- prot_df$organism_name[match(to_plot$organism_acc, prot_df$organism_acc)]
to_plot$organism_name <- paste0(to_plot$organism_name, " (", to_plot$organism_acc, ")")

ggplot(to_plot, aes(variable, organism_name, fill = variable, alpha = value)) + 
    geom_tile(colour = "gray50") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_minimal() + 
    labs(x = "", y = "") +
    coord_flip() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "None") +
    ggtitle("Absence Presence of pR1SE proteins", 
            subtitle = paste0("NCBI (added to 2 or more profiles n = ",nrow(to_plot), ")"))
    

ggsave("../../plots/NCBI_presence_absence_plot.png", plot = last_plot(), height = 8, width = 16)
fwrite(wide_df, "known_hosts_df.tsv")


# I need this for 03_retrieve_genomes.pu ----------------------------------
write(wide_df$organism_acc, file = "accs_to_download.txt")










# test --------------------------------------------------------------------

wide_df$first <- 0
wide_df$second <- 0

for(i in 1:nrow(wide_df)){
    wide_df$first[i] <- sum(wide_df$ORF6[i], wide_df$ORF8[i], wide_df$ORF10[i])
    wide_df$second[i] <- sum(wide_df$ORF17[i], wide_df$ORF21[i], wide_df$ORF23[i], wide_df$ORF24[i])
    
    
}







