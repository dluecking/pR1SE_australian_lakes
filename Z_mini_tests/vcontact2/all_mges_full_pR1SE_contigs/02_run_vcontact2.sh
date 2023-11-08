#!/bin/bash
#needs vContact2 env

vcontact2 --raw-proteins vcontact_input/protein_file.faa --rel-mode 'Diamond' --proteins-fp vcontact_input/gene_to_genome_df.csv --db 'ArchaeaViralRefSeq207-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin ~/bin/cluster_one-1.0.jar --output-dir vcontact2_out --verbose 
