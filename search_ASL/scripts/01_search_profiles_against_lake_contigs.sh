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
#SBATCH --array=1-12


source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate hmmer

echo ${SLURM_ARRAY_TASK_ID}

LAKE=`cat list_of_lakes.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

for ORF in `cat list_of_ORFs.txt`;
do
  if [ -f "../01_hmm_results/${LAKE}_${ORF}_domtblout.txt" ]; then
    ### Take action if $DIR exists ###
    echo "${ORF} was already search in ${LAKE}, doing nothing!"
  else
    ###  Control will jump here if $DIR does NOT exists
    hmmsearch \
      --domtblout  ../01_hmm_results/${LAKE}_${ORF}_domtblout.txt \
      --tblout ../01_hmm_results/${LAKE}_${ORF}.txt \
      --cpu 16 \
      ../input_profiles/${ORF}.faa.msa.hmm \
      ~/bioinf/australian_salt_lakes/tal/4_Dom/${LAKE}_metaspades.faa

  fi
done   
