#!/bin/bash

# needs vs2 env active for prodigal 
for file in ../contigs/*;
do
  contig=`basename $file`
  prodigal -p meta -i $file -f sco -o ../sco/${contig}.sco \
    -a ../genes/${contig}.faa
done

