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
conda activate iterative_hmm 

ORF=`cat ../../helper_files/list_of_ORFs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

python inflate_protein_to_protile_to_mafft.py \
  --protein_file ../05_curated_blastp_proteins/${ORF}_blastp_proteins.faa \
  --output_file ../06_hmm_iter1/${ORF}_iter1.msa \
  --database_path ~/bioinf/dbs/diamond_nr/nr.gz \
  --keep_intermediate

