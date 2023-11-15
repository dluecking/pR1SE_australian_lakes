# dlu 14.11.23
SER wants me to follow up on the self targeting found with the blast against the cripsr iphop database appraoch. For that I download the genomes of the apparently apHPV targeting chromosomes, via NCBI (search GCA_accession > download chromosome).
Then I predict the spacers by running the online CRIPSRCasFinder tool:
https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Executing

Then I select 'display spacers' and exclude evidence level 1 spacers.

Issue: out of 8 potential chromosomes, some are duplicate (they share the host), then out of the 6 unique, 2 are too fragmented (>100 contigs to be run), and 1 has no spacers (above level 1). This leaves 3 sets of spacers, which I will blast against apHPVs.

