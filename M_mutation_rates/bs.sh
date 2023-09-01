#!/bin/bash
echo "Note: This script is designed to run with the amount of memory detected by BBMap."
echo "      If Samtools crashes, please ensure you are running on the same platform as BBMap,"
echo "      or reduce Samtools' memory setting (the -m flag)."
samtools sort -m 10G -@ 3 bam/Gairdner_NODE_167.bam -o bam/Gairdner_NODE_167_sorted.bam
samtools index bam/Gairdner_NODE_167_sorted.bam
