#!/bin/bash

for ORF in `cat ../../helper_files/list_of_ORFs.txt`;
do
  python create_blacklist.py ../03_blastp_proteins/${ORF}_blastp_proteins.faa ../05_curated_blastp_proteins/${ORF}_blastp_proteins.faa >> ../accession_blacklist.txt
done



