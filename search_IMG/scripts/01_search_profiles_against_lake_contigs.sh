#!/bin/bash
# needs hmmer active
for ORF in `cat list_of_ORFs.txt`;
do
  if [ -f "../01_hmm_results/${ORF}_vs_IMG_domtblout.txt" ]; then
    ### Take action if $DIR exists ###
    echo "${ORF} was already searched, doing nothing!"
  else
    ###  Control will jump here if $DIR does NOT exists
    hmmsearch \
      --domtblout  ../01_hmm_results/${ORF}_vs_IMG_domtblout.txt \
      --tblout ../01_hmm_results/${ORF}_vs_IMG.txt \
      --cpu 16 \
      ../input_profiles/${ORF}.faa.msa.hmm \
      ~/bioinf/dbs/IMG_VR/IMG_VR_2022-12-19_7/IMGVR_all_proteins.faa.gz
  fi
done   
