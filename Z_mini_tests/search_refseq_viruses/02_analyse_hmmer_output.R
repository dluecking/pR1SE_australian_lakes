# author: dlueckin
# date: Tue May  9 13:14:32 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(tidyr)
library(dplyr)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# explore -----------------------------------------------------------------

big_df <- data.table()

for(file in list.files("hmm_out/", full.names = TRUE, pattern = "dom.out")){
    is_not_empty <- as.numeric(str_split(system(command = paste0("wc -l ", file), intern = TRUE), " ")[[1]][1]) > 13
    
    if(is_not_empty){
        domtbl <- fread(file, skip = 3, fill = TRUE) %>% 
            select(V1, V3, V4, V6, V7, V8, V16, V17) %>% 
            filter(V1 != "#")
        
        names(domtbl) <- c("hit", "hit_length", "query", "query_length", "evalue", "score", "alignment_start", "alignment_end")
        
        current_ORF <- str_extract(file, pattern = "ORF[^\\_]*")
        
        # calc alignment length
        domtbl$aln_length <- domtbl$alignment_end - domtbl$alignment_start
        
        # prepare data for later filtering
        domtbl$query_length <- as.numeric(domtbl$query_length)
        domtbl$score <- as.numeric(domtbl$score)
        
        domtbl <- domtbl

        big_df <- rbind(big_df, domtbl)
    }
}

rm(domtbl, current_ORF)

# add uvig
big_df$uvig <- str_extract(big_df$hit, pattern = "^[^|]*")
# add host taxonomy
all_host_info <- fread("../../../../bioinf/dbs/IMG_VR/IMG_VR_2022-12-19_7/IMGVR_all_Host_information.tsv")
big_df$host_taxonomy <- all_host_info$`Host taxonomy prediction`[match(big_df$uvig, all_host_info$UVIG)]

# add phage taxonomy
all_seq_info <- fread("../../../../bioinf/dbs/IMG_VR/IMG_VR_2022-12-19_7/IMGVR_all_Sequence_information.tsv")
big_df$virus_taxonomy <- all_seq_info$`Taxonomic classification`[match(big_df$uvig, all_seq_info$UVIG)]
# add virus completeness
big_df$virus_completeness <- all_seq_info$`Estimated completeness`[match(big_df$uvig, all_seq_info$UVIG)]
big_df$virus_quality <- all_seq_info$`MIUViG quality`[match(big_df$uvig, all_seq_info$UVIG)]

fwrite(big_df, "hit_taxonomy.csv")
