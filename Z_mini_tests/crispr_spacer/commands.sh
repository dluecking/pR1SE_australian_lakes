for i in contigs/*; 
do
    blastn -query $i \
    -db ~/bioinf/dbs/iphop_db/Test_db/db/All_CRISPR_spacers_nr_clean \
    -out ${i}_vs_test_db.out -outfmt 6; 
done