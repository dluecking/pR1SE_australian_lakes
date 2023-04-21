#!/bin/sh
# properties = {"type": "single", "rule": "iterative_blastp", "local": false, "input": ["initial_proteins/ORF6_initial_proteins.faa"], "output": ["blastp_proteins/ORF6_blastp_proteins.faa"], "wildcards": {"ORF": "ORF6"}, "params": {"db": "/bioinf/home/dlueckin/dbs/diamond_nr", "evalue": 1e-05, "bitscore": 50}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp", "mem": "50GB", "threads": 16, "time": "10:00:00"}, "jobid": 0, "cluster": {}}
 cd /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters && \
/home/dlueckin/bin/miniconda3/envs/sm/bin/python3.9 \
-m snakemake blastp_proteins/ORF6_blastp_proteins.faa --snakefile /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/Snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/.snakemake/tmp.yy2s2ztz' 'initial_proteins/ORF6_initial_proteins.faa' --latency-wait 50 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules iterative_blastp --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /home/dlueckin/bin/miniconda3/envs/sm/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/.snakemake/tmp.yy2s2ztz/0.jobfinished || (touch /home/dlueckin/projects/pR1SE_australian_lakes/generate_protein_clusters/.snakemake/tmp.yy2s2ztz/0.jobfailed; exit 1)

