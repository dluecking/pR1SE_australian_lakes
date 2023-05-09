#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=hmmer   # Job name
#SBATCH --output=out/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=1-00:00:00   # walltime
#SBATCH --mem=100G   # memory per NODE
#SBATCH --ntasks=16
#SBATCH --partition=CLUSTER
#SBATCH --array=2-8


source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate hmmer

echo ${SLURM_ARRAY_TASK_ID}

ORF=`cat list_of_ORFs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

if [ -f "../08_hmmer_results/${ORF}.txt" ]; then
  ### Take action if $DIR exists ###
  echo "${ORF} already exists, doing nothing!"
else
  ###  Control will jump here if $DIR does NOT exists
  hmmsearch \
    --domtblout ../08_hmmer_results/${ORF}_domtblout.txt \
    --tblout ../08_hmmer_results/${ORF}.txt \
    --cpu 16 \
    ../07_curated_profiles/${ORF}_blastp_proteins_curated.faa.msa.hmm \
    ~/bioinf/dbs/diamond_nr/nr.gz
fi
 
