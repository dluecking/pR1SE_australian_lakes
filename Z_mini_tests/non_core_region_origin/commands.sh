# make blast db
makeblastdb -in combined_all_contigs.fasta -dbtype nucl -title combined_all_contigs -out combined_all_contigs

# blast contigs against db
for i in ../../../FINAL_complete_relatives/*; do blastn -query $i -db ../all_halobacteria/blast_db/combined_all_contigs -out `basename $i`_vs_all_halos.out -outfmt 6; done
