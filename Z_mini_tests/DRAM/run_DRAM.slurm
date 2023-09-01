#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=DRAM_pR1SE   # Job name
#SBATCH --output=out/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=4-00:00:00   # walltime
#SBATCH --mem=200G   # memory per NODE
#SBATCH --ntasks=16
#SBATCH --partition=CLUSTER    # Partition
#SBATCH --array=2-38 # max 38


source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate DRAM

echo ${SLURM_ARRAY_TASK_ID}
acc=`cat list_of_mags.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

DRAM.py annotate -i ../../FINAL_complete_relatives/${acc}.fasta -o DRAM_out/${acc} --threads 16
