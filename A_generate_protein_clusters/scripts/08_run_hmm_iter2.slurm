#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=generate_protein_cluster_iter2   # Job name
#SBATCH --output=out/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=2-00:00:00   # walltime
#SBATCH --mem=100G   #memory per NODE
#SBATCH --cpus-per-task=16
#SBATCH --partition=CLUSTER
#SBATCH --array=1-7

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate iterative_hmm

ORF=$(cat "../../helper_files/list_of_ORFs.txt" | sed -n "${SLURM_ARRAY_TASK_ID}p")

python inflate_protein_to_protile_to_mafft.py \
	--protein_file "../07_curated_iter1_proteins/${ORF}_proteins.fasta" \
	--output_file "../08_hmm_iter2/${ORF}_iter2.msa" \
	--database_path ~/bioinf/dbs/diamond_nr/nr.gz \
	--keep_intermediate
