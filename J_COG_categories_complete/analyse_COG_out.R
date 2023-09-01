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

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# load important dfs ------------------------------------------------------

manual_region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                             sheet = "pR1SE_relatives"))
manual_region_df <- manual_region_df %>% 
    filter(contig_category == "complete") %>% 
    filter(!is.na(manual_min_orf))



# load prot df ------------------------------------------------------------

prot_df <- fread("../G_annotate_complete_relatives/scripts/improved_gene_df_with_DRAM_and_nr.tsv")


# load blast_df -----------------------------------------------------------

blast_df <- rbindlist(lapply(list.files("blast_out", full.names = TRUE), fread))


# add COG db, based on hit ------------------------------------------------

def <- fread("../../misc/test_new_COG/cog-20.def.tab")
# cog needs to be treated with sed -i 's/,,,,,/,,,,/' cog-20.cog.csv
cog <- fread("../../misc/test_new_COG/cog-20.cog.csv", fill = T, sep = ",") 

blast_df$COG <- cog$V7[match(blast_df$V2, cog$V3)]


# add actual letter per hit -----------------------------------------------

blast_df$letter <- def$V2[match(blast_df$COG, def$V1)]


# add orf number, extract contig ------------------------------------------

blast_df$orf <- as.numeric(str_extract(blast_df$V1, pattern = "\\d*$"))
blast_df$contig <- str_remove(blast_df$V1, pattern = "\\_\\d*$")


# filter blast_df by 10-5  ------------------------------------------------

blast_df <- blast_df %>% 
    filter(V11 <= 10^-5)


# select relevant columns -------------------------------------------------

blast_df <- blast_df %>% 
    select(V1, COG, letter, orf, contig)
names(blast_df) <- c("id", "COG", "letter", "orf", "contig")


# add prot_df annotation --------------------------------------------------

blast_df$annotation <- prot_df$annot[match(blast_df$id, prot_df$id)]


# filter out all non-relevant contigs -------------------------------------

blast_df <- blast_df %>% 
    filter(contig %in% manual_region_df$contig)


# label each gene either "core", "inside", "outside" ----------------------

# everything gets unknown
blast_df$type <- "unknown"

for(i in 1:nrow(blast_df)){
    # current contig
    CONTIG <- blast_df$contig[i]
    ORF <- blast_df$orf[i]
    CIRCULARITY <- manual_region_df$genome_segment[manual_region_df$contig == CONTIG]
    STATE <- manual_region_df$integrated_state[manual_region_df$contig == CONTIG]
    
    # if we are within ORF10 and ORF 17
    CURRENT_ORF10_pos <- prot_df %>% 
        filter(contig == CONTIG, annot == "ORF10.faa") %>% 
        select(`orf#`) %>% 
        pull()
    
    CURRENT_ORF17_pos <- prot_df %>% 
        filter(contig == CONTIG, annot == "ORF17.faa") %>% 
        select(`orf#`) %>% 
        pull()
    # if we have an error
    if(length(CURRENT_ORF17_pos) > 1 | length(CURRENT_ORF10_pos) > 1){
        
        if(CONTIG == "LWLN01000003.1"){
            #... its because we are LWLN01000003
            # then take the first region
            CURRENT_ORF10_pos <- 16
            CURRENT_ORF17_pos <- 21
        }
        
        if(CONTIG == "CP100004.1"){
            # or we are in CP100004.1
            # then we just give back an error, too complex to fix
            blast_df$type[i] <- "error"
            next
        }
        
        
    }
    
    # if we are between these two positions, we are inside
    if(between(ORF, left = CURRENT_ORF10_pos + 1, right = CURRENT_ORF17_pos - 1)){
        blast_df$type[i] <- "inside"
        next
    }
    
    # if we are inside region borders -4 and +4 
    # (to re-adjust the region border finding we did)
    # we get "core"
    REGION_MIN <- manual_region_df$manual_min_orf[manual_region_df$contig == CONTIG] + 4
    REGION_MAX <- manual_region_df$manual_max_orf[manual_region_df$contig == CONTIG] - 4
    
    if(between(ORF, left = REGION_MIN, REGION_MAX)){
        # we cant hit inside, otherwise we would have left the loop at last "next"
        blast_df$type[i] <- "core" 
        next
    }
    
    # last check: if we are outside, we need to be on a "cirucular" contig that is not a main chrom, otherwise
    if(STATE != "main_chromosome"){
        #CIRCULARITY == "circular" &
        # we would be outside the region
        blast_df$type[i] <- "outside" 
        next
    }
    
}


# the following line tells me, that we cant analyse inside and core, since not enough data
# we can do "outside", if we include non-circular plasmids
table(blast_df$type)


# prepare for plotting ----------------------------------------------------

blast_df_filtered <- blast_df %>% 
    filter(type == "outside")

template_df <- fread("../../../template_classifier_count.tsv")

for(i in 1:nrow(template_df)){
    current_letter <- template_df$LETTER[i]
    template_df$COUNT[i] <- blast_df_filtered %>% 
        filter(str_detect(string = letter, pattern = current_letter)) %>% 
        nrow() %>% 
        unlist()
}
template_df$perc <- template_df$COUNT / sum(template_df$COUNT) * 100


# plot
ggplot(template_df, aes(x = SHORT_DESCRIPTION, y = perc)) +
    geom_bar(stat = 'identity', color = "black", fill = "grey", alpha = 0.9) +
    theme_minimal() +
    ggtitle("COG category of genes OUTSIDE",
            subtitle = "outside: plasmids only, outside of pR1SE_regions") +
    theme(axis.text.x = element_text(angle = 45 ,  hjust=1),
          plot.margin=unit(c(0.5,0.5,0.5,l = 2), "cm"))
ggsave(plot = last_plot(), filename = "../plots/COG_outside_pR1SE_regions.pdf", height = 4, width = 6)
