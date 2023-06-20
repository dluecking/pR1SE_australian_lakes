## How did I get to the minscores for each ORF?
1. combine ALL blastp-found proteins 
`cat 03_blastp_proteins/*.faa$ > combined_input_proteins.faa`

2. create an hmm for each set of curated blastp-found proteins
```
cd ../05_curated_blastp_proteins
for i in *; do mafft --localpair --reorder $i > ../00_determine_minscore/curated_blastp_msa/${i}.msa; done

cd ../00_determine_minscore/curated_blastp_msa/
for i in *; do hmmbuild --amino ../curated_blastp_hmm/${i}.hmm $i; done
```
3. search the profiles against the blastp-found uncurated proteins
`for i in *; do hmmsearch --domtblout ../hmm_out/${i}.domtblout --cpu 40 $i ../all_blastp_proteins.faa ; done`
