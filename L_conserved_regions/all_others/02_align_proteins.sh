# needs mafft active
for i in proteins/*.faa
do 
    acc=`basename $i`
    mafft --localpair --reorder $i > alignments/${acc}.msa
done
