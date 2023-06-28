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
#SBATCH --mem=40G   # memory per NODE
#SBATCH --ntasks=22
#SBATCH --partition=SMALLMEM    # Partition


source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate iphop_mamba
ls /bioinf/transfer/tal/iphop_db/db/
iphop predict --fa_file ../FINAL_complete_relatives/CP031299.1.fasta --db_dir /bioinf/transfer/tal/iphop_db/db/ --out_dir here
