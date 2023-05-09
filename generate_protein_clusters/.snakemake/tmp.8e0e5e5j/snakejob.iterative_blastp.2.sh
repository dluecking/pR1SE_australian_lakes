#!/bin/sh
# properties = {"type": "single", "rule": "iterative_blastp", "local": false, "input": ["initial_proteins/ORF8_initial_proteins.faa"], "output": ["blastp_proteins/ORF8_blastp_proteins.faa"], "wildcards": {"ORF": "ORF8"}, "params": {"db": "/bioinf/home/dlueckin/dbs/diamond_nr/nr.dmnd", "evalue": 1e-05, "bitscore": 50}, "log": ["log/ORF8_test.log"], "threads": 1, "resources": {"tmpdir": "/tmp", "mem": "50GB", "threads": 16, "time": "10:00:00", "partition": "CLUSTER"}, "jobid": 2, "cluster": {}}
 cd /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters && \
/home/dlueckin/bin/miniconda3/envs/sm/bin/python3.9 \
-m snakemake blastp_proteins/ORF8_blastp_proteins.faa --snakefile /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/Snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/.snakemake/tmp.8e0e5e5j' 'initial_proteins/ORF8_initial_proteins.faa' --latency-wait 50 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules iterative_blastp --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /home/dlueckin/bin/miniconda3/envs/sm/bin \
--mode 2  --use-conda  --conda-frontend mamba  --conda-base-path /home/dlueckin/bin/miniconda3  --default-resources "tmpdir=system_tmpdir"  && touch /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/.snakemake/tmp.8e0e5e5j/2.jobfinished || (touch /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/.snakemake/tmp.8e0e5e5j/2.jobfailed; exit 1)

