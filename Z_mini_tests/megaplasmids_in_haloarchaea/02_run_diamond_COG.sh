#!/bin/bash

# needs diamond env active
for i in genes/*;
do
    acc=`basename $i`
    diamond blastp \
    --query $i \
    --db ~/bioinf/dbs/COGdb/dmnd_db/combined_COG.dmnd  \
    --out blast_out/${acc}.out -f 6 --max-target-seqs 1
done
