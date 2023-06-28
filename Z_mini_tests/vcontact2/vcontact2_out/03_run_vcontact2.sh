#!/bin/bash
#needs vContact2 env

vcontact2 --raw-proteins ../input_files/protein_file.faa --rel-mode 'Diamond' --proteins-fp ../input_files/gene_to_genome_df.csv --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin ~/bin/cluster_one-1.0.jar --output-dir ../vcontact2_out
