#!/bin/bash
# needs DRAM

for i in DRAM_out/*;
do
    acc=$(basename $i)
    echo $acc

    TRNA_FILE=DRAM_out/${acc}/trnas.tsv

    if [ ! -f $TRNA_FILE ]; then
        DRAM.py distill \
            -i DRAM_out/${acc}/annotations.tsv \
            -o distill_out/${acc} 
    else
        DRAM.py distill \
            -i DRAM_out/${acc}/annotations.tsv \
            -o distill_out/${acc} \
            --trna_path $TRNA_FILE 
    fi
done