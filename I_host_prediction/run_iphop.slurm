#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=iphop   # Job name
#SBATCH --output=out/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=4-00:00:00   # walltime
#SBATCH --mem=200G   # memory per NODE
#SBATCH --ntasks=16
#SBATCH --partition=SMALLMEM    # Partition
#SBATCH --array=1-18

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate iphop_mamba

echo ${SLURM_ARRAY_TASK_ID}
acc=`cat new_list_of_hosts_to_predict.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

iphop predict \
--fa_file ../FINAL_complete_relatives/${acc}.fasta \
--db_dir /bioinf/transfer/tal/iphop_db/ \
--out_dir iphop_out/${acc}
