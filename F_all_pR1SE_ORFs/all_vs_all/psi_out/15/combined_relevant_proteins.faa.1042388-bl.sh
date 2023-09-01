#!/bin/sh
#PBS -v PATH
#$ -v PATH


para=$1
cd /home/domi/Documents/psi_out
./combined_relevant_proteins.faa.1042388-bl.pl 0 $para &
wait

