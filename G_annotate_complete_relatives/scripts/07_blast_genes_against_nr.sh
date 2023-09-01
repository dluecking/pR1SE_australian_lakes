#!/bin/bash
# needs diamond active

diamond blastp \
--query ../all_pR1SE_relatives_genes.faa \
--db ~/bioinf/dbs/diamond_nr/nr.dmnd \
-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
--out ../diamond_blastp_nr_out/all_pR1SE_relatives_genes_vs_nr.out \
--max-target-seqs 3 \
--evalue 0.00001 \
--memory-limit 50G \
--threads 40 \
--fast
