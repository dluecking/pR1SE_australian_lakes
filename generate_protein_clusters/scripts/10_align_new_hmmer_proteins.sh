#!/bin/bash
cd ../10_hmmer_proteins/
for i in *; do mafft --localpair --reorder $i > ../11_hmmer_alignments/${i}.msa ; done
