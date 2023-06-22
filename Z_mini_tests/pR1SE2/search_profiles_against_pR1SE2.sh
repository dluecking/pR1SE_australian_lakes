#!/bin/bash -l

for ORF in `cat ../../helper_files/list_of_ORFs.txt`;
do
  if [ -f "../01_hmm_results/pR1SE2_${ORF}_domtblout.txt" ]; then
    ### Take action if $DIR exists ###
    echo "${ORF} was already search in pR1SE2, doing nothing!"
  else
    ###  Control will jump here if $DIR does NOT exists
    hmmsearch \
      --domtblout  01_hmm_results/pR1SE2_${ORF}_domtblout.txt \
      --tblout 01_hmm_results/pR1SE2_${ORF}.txt \
      --cpu 16 \
      ../../A_generate_protein_clusters/FINAL_profiles/${ORF}.faa.msa.hmm \
      pR1SE2_proteins.faa
  fi
done   
