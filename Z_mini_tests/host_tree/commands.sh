# the gtdbk command from seba
gtdbtk classify_wf --cpus 24 --genome_dir host_genomes/ --out_dir gtdbtk_out --extension fn

# did it online but they said they ran this command
path_to_iqtree -s gtdbtk.ar53.user_msa.fasta.gz -m TEST -bb 1000 -alrt 1000
