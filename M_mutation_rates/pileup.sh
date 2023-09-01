#/bin/bash
# needs mapping active

samtools mpileup \
    -f contigs/Eyre_NODE_236_length_33155_cov_19.932.fasta \
    bam/Eyre_NODE_236_sorted_mappedonly.bam \
    > pileup/Eyre_NODE_236_sorted_mappedonly.bam.pileup

samtools mpileup \
    -f contigs/Eyresouth_NODE_101_length_48437_cov_7.17035.fasta \
    bam/Eyresouth_NODE_101_sorted_mappedonly.bam \
    > pileup/Eyresouth_NODE_101_sorted_mappedonly.bam.pileup

samtools mpileup \
    -f contigs/Frome_NODE_383_length_24936_cov_4.316346.fasta \
    bam/Frome_NODE_383_sorted_mappedonly.bam \
    > pileup/Frome_NODE_383_sorted_mappedonly.bam.pileup

samtools mpileup \
    -f contigs/Gairdner_NODE_107_length_49922_cov_74.384603.fasta \
    bam/Gairdner_NODE_107_sorted_mappedonly.bam \
    > pileup/Gairdner_NODE_107_sorted_mappedonly.bam.pileup

samtools mpileup \
    -f contigs/Gairdner_NODE_167_length_41210_cov_9.888325.fasta \
    bam/Gairdner_NODE_167_sorted_mappedonly.bam \
    > pileup/Gairdner_NODE_167_sorted_mappedonly.bam.pileup

