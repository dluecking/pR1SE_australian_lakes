# needs hmmer env active
for ORF in `cat ../../helper_files/list_of_ORFs.txt`;
do
  for gene_file in ../genes/*;
  do
    acc=`basename $gene_file`
    if [ -f "../hmm_results/${acc}_vs_${ORF}_domtblout.txt" ]; then
      ### Take action if $DIR exists ###
      echo "${ORF} was already search in ${acc}, doing nothing!"
    else
      ###  Control will jump here if $DIR does NOT exists
      hmmsearch \
        --domtblout  "../hmm_results/${acc}_vs_${ORF}_domtblout.txt" \
        --cpu 16 \
        ../../A_generate_protein_clusters/FINAL_profiles/${ORF}.faa.msa.hmm \
        $gene_file
    fi
  done
done

