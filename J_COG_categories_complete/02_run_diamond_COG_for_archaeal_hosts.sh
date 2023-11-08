#!/bin/bash

# needs diamond env active
diamond blastp \
    --query archaeal_proteins/all_halo_proteins_combined.faa \
    --db ~/bioinf/dbs/COGdb/dmnd_db/combined_COG.dmnd  \
    --out blast_out_archaeal_proteins/blast_vs_COG.out -f 6 --max-target-seqs 1 \
    --threads 24
