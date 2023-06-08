for i in interpro_in/*
do 
  name=`basename $i`
  python ~/bin/webservice-clients/python/iprscan5.py \
    --sequence=${i} \
    --outfile interpro_out/${name}.out \
    --email=dlueckin@mpi-bremen.de 
done
