#!/bin/bash
#needs vContact2 env

vcontact2 --raw-proteins protein_file.faa --rel-mode 'Diamond' --proteins-fp gene_to_genome_df.csv --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin ~/bin/cluster_one-1.0.jar --output-dir vcontact2_out --verbose 
