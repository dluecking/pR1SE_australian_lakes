#!/bin/bash
# needs hmmer env active
for i in ../13_curated_alignments/*; 
do 
  hmmbuild --amino ../14_curated_profiles/${i}.hmm $i; 
done

