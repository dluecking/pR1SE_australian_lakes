# needs conda tRNAscanSE active

for i in ../../FINAL_complete_relatives/*; do acc=`basename ${i}`; tRNAscan-SE -A ${i} -o ${acc}.out -m ${acc}.summary; done
