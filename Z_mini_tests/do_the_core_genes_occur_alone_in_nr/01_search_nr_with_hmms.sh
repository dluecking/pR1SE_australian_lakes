#!/bin/bash
# needs iterative_hmm

for i in ../../A_generate_protein_clusters/FINAL_profiles/*; 
do 
    ORF=`basename $i`;
    hmmsearch --domtblout ${ORF}_vs_nr.domtblout --cpu 24 $i ~/bioinf/dbs/diamond_nr/nr.gz; 
done



