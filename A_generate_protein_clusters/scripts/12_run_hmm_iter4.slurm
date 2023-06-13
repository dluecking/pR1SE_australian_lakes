#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=generate_protein_cluster_iter3   # Job name
#SBATCH --output=out/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=2-00:00:00   # walltime
#SBATCH --mem=100G   #memory per NODE
#SBATCH --cpus-per-task=16
#SBATCH --partition=CLUSTER
#SBATCH --array=1-2 # only ORFS 8 and 10 are to be processed, the rest is done.

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate iterative_hmm

ORF=$(cat "../../helper_files/list_of_ORFs.txt" | sed -n "${SLURM_ARRAY_TASK_ID}p")

python inflate_protein_to_protile_to_mafft.py \
	--protein_file "../11_curated_iter3_proteins/${ORF}_proteins.fasta" \
	--output_file "../12_hmm_iter4/${ORF}_iter3.msa" \
	--database_path ~/bioinf/dbs/diamond_nr/nr.gz \
	--keep_intermediate
