# author: dlueckin
# date: Mon Sep 18 10:15:38 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gggenes)
library(googlesheets4)


# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# read google sheet -------------------------------------------------------
CONTIGS_TO_PLOT <- c("LWLN01000003.1", "Gairdner_NODE_107_length_49922_cov_74.384603", 
                     "NC_013748.1", "SJER01000008.1", "KX687704.1")

region_df <- as.data.table(read_sheet("https://docs.google.com/spreadsheets/d/1J9eVHjpSJsBOVN_gavMRTGaNeT7iXZSdowyOkn9GUsA/edit#gid=1238763630",
                                      sheet = "pR1SE_relatives"))

region_df <- region_df %>% 
    filter(contig %in% CONTIGS_TO_PLOT)


# generate gene_df --------------------------------------------------------

gene_df <- data.table()

for(contig in CONTIGS_TO_PLOT){
    file <- paste0("../E_check_synteny/sco/", contig, ".fasta.sco")
    df <- fread(file, skip = 2, sep = "_")
    names(df) <- c("orf#", "start", "end", "strand")
    df$strand[df$strand == "+"] <- TRUE
    df$strand[df$strand == "-"] <- FALSE
    df$`orf#` <- as.numeric(str_remove(df$`orf#`, ">"))
    df$contig <- contig
    
    gene_df <- rbind(gene_df, df)
}
rm(contig, file, df)

gene_df$id <- paste0(gene_df$contig, "_", gene_df$`orf#`)



# add correct annotation --------------------------------------------------

annot_df <- fread("../E_check_synteny/scripts/cleaned_gene_df_complete.tsv")

gene_df$annot <- annot_df$annot[match(gene_df$id, annot_df$id)]
gene_df$annot[is.na(gene_df$annot)] <- NA


# only keep genes within region -------------------------------------------


gene_df_filtered <- data.table()

for(i in 1:nrow(region_df)){
    
    filtered <- gene_df %>% 
        filter(contig == region_df$contig[i]) %>% 
        filter(between(x = `orf#`, left = region_df$manual_min_orf[i], right = region_df$manual_max_orf[i]))
    
    gene_df_filtered <- rbind(gene_df_filtered, filtered)
}
gene_df <- gene_df_filtered
rm(gene_df_filtered, i, filtered)


# cosmetics ---------------------------------------------------------------

# shorter label
gene_df$contig <- str_remove(gene_df$contig, "\\_length.*")

# add the two other ORFs
gene_df$annot[which(gene_df$annot == "ORF6") + 1] <- "ORF7"
gene_df$annot[which(gene_df$annot == "ORF24") + 1] <- "ORF25"




# plot --------------------------------------------------------------------


ggplot(gene_df, aes(xmin = start, xmax = end, y = contig, fill = annot, label = annot)) +
    geom_gene_arrow() +
    geom_gene_label(align = "left") +
    facet_wrap(~factor(contig, levels = c("SJER01000008.1", 
                                          "LWLN01000003.1", 
                                          "KX687704.1",
                                          "NC_013748.1",  
                                          "Gairdner_NODE_107_length_49922_cov_74.384603")),
               scales = "free", ncol = 1, ) +
    theme_genes() +
    scale_fill_manual(values = c("ORF6" = "#f94144",
                                 "ORF8" = "#f3722c",
                                 "ORF10" = "#f8961e",
                                 "ORF17" = "#f9c74f",
                                 "ORF21" = "#90be6d",
                                 "ORF23" = "#43aa8b",
                                 "ORF24" = "#577590",
                                 "unknown" = "#FBFEF9",
                                 "ORF7" = "#D0D3D6",
                                 "ORF25" = "#D0D3D6")) +
    theme(legend.position = "None") +
    ylab("") +
    xlab("Position [bp]") +
    theme(
        # panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ) +
    theme(plot.margin=unit(c(0.5,0.5,r = 0.5, 0.1), "cm"))

# ggsave(plot = last_plot(), file = "genome_maps_5_representatives.png", height = 5, width = 8)
# ggsave(plot = last_plot(), file = "genome_maps_5_representatives.svg", height = 5, width = 8)



# add transparency values -------------------------------------------------

blast_df <- rbindlist(lapply(list.files("../F_all_pR1SE_ORFs/all_vs_all/blast_out/", full.names = TRUE), fread))
names(blast_df) <- unlist(str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", " "))

relevant_ids <- gene_df %>% 
    filter(!is.na(annot)) %>% 
    select(id) %>% 
    pull()

relevant_blast_df <- blast_df %>% 
    filter(qseqid %in% relevant_ids & sseqid %in% relevant_ids)
relevant_blast_df <- relevant_blast_df %>% 
    select(qseqid, sseqid, pident)

relevant_blast_df$q_annot <- gene_df$annot[match(relevant_blast_df$qseqid, gene_df$id)]
relevant_blast_df$s_annot <- gene_df$annot[match(relevant_blast_df$sseqid, gene_df$id)]

relevant_blast_df <- relevant_blast_df %>% 
    filter(q_annot != "unknown") %>% 
    filter(s_annot != "unknown") %>% 
    filter(sseqid != qseqid)
