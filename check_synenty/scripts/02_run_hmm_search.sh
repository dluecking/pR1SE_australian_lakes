
for ORF in `cat list_of_ORFs.txt`;
do
  for gene_file in ../genes/complete/*;
  do
    acc=`basename $gene_file`
    if [ -f "../hmm_results/complete/${acc}_vs_${ORF}_domtblout.txt" ]; then
      ### Take action if $DIR exists ###
      echo "${ORF} was already search in ${acc}, doing nothing!"
    else
      ###  Control will jump here if $DIR does NOT exists
      hmmsearch \
        --domtblout  "../hmm_results/complete/${acc}_vs_${ORF}_domtblout.txt" \
        --cpu 16 \
        ../input_profiles/${ORF}.faa.msa.hmm \
        $gene_file
    fi
  done
done

