#!/bin/bash

# needs DRAM env active
for file in ../../FINAL_complete_relatives/*;
do
    acc=`basename $file`
    DRAM.py annotate \
    -i $file \
    -o ../DRAM_out/${acc} \
    --min_contig_size 10000 \
    --prodigal_mode meta


    DRAM.py distill \
    -i ../DRAM_out/${acc}/annotations.tsv \
    -o ../DRAM_distilled/${acc} \
    --trna_path ../DRAM_out/${acc}/trnas.tsv \
    --rrna_path ../DRAM_out/${acc}/rrnas.tsv

done