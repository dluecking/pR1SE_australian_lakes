#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=generate_protein_cluster   # Job name
#SBATCH --output=out/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=2-00:00:00   # walltime
#SBATCH --mem=100G   #memory per NODE
#SBATCH --cpus-per-task=24
#SBATCH --partition=CLUSTER
#SBATCH --array=1-7

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate diamond

echo ${SLURM_ARRAY_TASK_ID}

ORF=`cat ../../helper_files/list_of_ORFs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

if [ -f ../03_blastp_proteins/${ORF}_blastp_proteins.faa ];
then
  echo "${ORF} already processed, doing nothing."
else
  echo "running: $ORF"
  python iterative_blast.py ../01_input_proteins/${ORF}_initial_proteins.faa \
    ../03_blastp_proteins/${ORF}_blastp_proteins.faa \
    --db ~/bioinf/dbs/diamond_nr/nr.dmnd \
    --evalue 0.00001 --bitscore 50 --threads 24
fi


