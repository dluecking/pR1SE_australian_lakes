## pR1SE in Australian Salt Lakes
### Task Description Susanne Email

We will use a core set of conserved proteins for the detection of new pR1SE relatives in the Australian salt lakes and public databases?

If we identified like 30-50 (?) possible candidates, we do the follwing:
1. Are there more possible "core proteins", if yes, which ones?
2. Are these core proteins on a main chromosome or on a plasmid?
3. What replication protein are found?
4. Who is the host?
5. Phylogeny of the relatives.

### Method description
1. generate a trusted database for each conserved proteins
    1. query public databases / Australian salt lakes with the proteins and find potentially more (evalue 10^-5 and score > 50).
    2. repeat this iteratively
    3. create a HMM (MAFFT for alignment and HMMER for hmm creation)
    4. again, iteratively search ASL and public databases thus increasing the number of proteins in a cluster even more (cutoffs are dependend on a search hmm-vs-pfam)
2. then screen binned/not-binned contigs and select ones with hits >4000 bp and at least 4 (might be more for us, since for 4/4 for pleos, so we have to go higher?)
3. further downstream quality control 

### Notes from Tomas
- the concept is: lose addition of homologs then later on removal of short proteins and more distant proteins (see step 5, where we set the minimum score to be included equal to the maximum score of false positives + 10)
- he has some pre-sorted assemblies where the contigs are sorted into plasmid / 


