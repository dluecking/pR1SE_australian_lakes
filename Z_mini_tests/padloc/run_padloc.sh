# needs padloc env


for i in `cat ../../helper_files/list_of_apHPVs.txt`;
do
	padloc --faa ../../FINAL_complete_relatives_proteins/aa/${i}.faa \
	--gff ../../FINAL_complete_relatives_proteins/gff/${i}.gff \
	--fix-prodigal --cpu 40 --outdir out
done
