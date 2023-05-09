#!/bin/bash
# needs hmmer env active
for i in *; do hmmbuild --amino ../07_curated_profiles/${i}.hmm $i; done

