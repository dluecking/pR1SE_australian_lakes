for i in proteins/*
do 
    acc=`basename $i`
    diamond blastp -q proteins/${acc} --db diamond_db/combined_relevant_proteins.faa.dmnd --out  blast_out/${acc}.out
done

