#!/bin/bash
cd ../12_curated_proteins/
for i in *; do mafft --localpair --reorder $i > ../13_curated_alignments/${i}.msa ; done
