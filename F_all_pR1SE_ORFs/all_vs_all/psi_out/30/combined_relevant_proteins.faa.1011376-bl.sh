#!/bin/sh
#PBS -v PATH
#$ -v PATH


para=$1
cd /home/domi/Documents
./combined_relevant_proteins.faa.1011376-bl.pl 0 $para &
wait

