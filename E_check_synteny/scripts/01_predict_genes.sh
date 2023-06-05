#!/bin/bash

# needs vs2 env active for prodigal 
for file in ../contigs/complete/*;
do
  contig=`basename $file`
  prodigal -i $file -f sco -o ../sco/complete/${contig}.sco \
    -a ../genes/complete/${contig}.faa
done


