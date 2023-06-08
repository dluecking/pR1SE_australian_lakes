#!/bin/bash
# needs mafft env active
cd ../03_blastp_proteins/
for i in *.faa; do mafft --localpair --reorder $i > ../04_blastp_alignments/${i}.msa ; done
