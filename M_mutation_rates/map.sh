#!/bin/bash
# needs bbmap active

bbmap.sh \
  in=reads/eyre_R1_trimmed_renamed.fq.gz \
  in2=reads/eyre_R2_trimmed_renamed.fq.gz \
  ref=contigs/Eyre_NODE_236_length_33155_cov_19.932.fasta \
  out=bam/Eyre_NODE_236.bam \
  slow k=12 bamscript=bs.sh; sh bs.sh nodisk

bbmap.sh \
  in=reads/eyresouth_R1_trimmed_renamed.fq.gz \
  in2=reads/eyresouth_R2_trimmed_renamed.fq.gz \
  ref=contigs/Eyresouth_NODE_101_length_48437_cov_7.17035.fasta \
  out=bam/Eyresouth_NODE_101.bam \
  slow k=12 bamscript=bs.sh nodisk; sh bs.sh

bbmap.sh \
  in=reads/frome_R1_trimmed_renamed.fq.gz \
  in2=reads/frome_R2_trimmed_renamed.fq.gz \
  ref=contigs/Frome_NODE_383_length_24936_cov_4.316346.fasta \
  out=bam/Frome_NODE_383.bam \
  slow k=12 bamscript=bs.sh nodisk; sh bs.sh

bbmap.sh \
  in=reads/gairdner_R1_trimmed_renamed.fq.gz \
  in2=reads/gairdner_R2_trimmed_renamed.fq.gz \
  ref=contigs/Gairdner_NODE_107_length_49922_cov_74.384603.fasta \
  out=bam/Gairdner_NODE_107.bam \
  slow k=12 bamscript=bs.sh nodisk; sh bs.sh 

bbmap.sh \
  in=reads/gairdner_R1_trimmed_renamed.fq.gz \
  in2=reads/gairdner_R2_trimmed_renamed.fq.gz \
  ref=contigs/Gairdner_NODE_167_length_41210_cov_9.888325.fasta \
  out=bam/Gairdner_NODE_167.bam \
  slow k=12 bamscript=bs.sh nodisk; sh bs.sh
