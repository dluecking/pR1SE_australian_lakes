#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=run_interpro   # Job name
#SBATCH --output=out/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=2-00:00:00   # walltime
#SBATCH --mem=10G   #memory per NODE
#SBATCH --cpus-per-task=4
#SBATCH --partition=CLUSTER
#SBATCH --array=1-38 # 38

FILE=$(cat "list_of_gene_files.txt" | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -f "../interpro_out/${FILE}.out" ]; then
	echo "${FILE} already exists, doing nothing"
else
	~/bioinf/dbs/interproscan-5.62-94.0/interproscan.sh \
	-o ../interpro_out/${FILE}.out \
	-cpu 24 \
	-f tsv \
	-i ../genes/${FILE} \
	--tempdir ../interpro_tmp/${FILE}
fi
