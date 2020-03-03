#!/bin/bash

set -e

bp=76

# We assume the following sharing matrix:
#   (A,B): 2
#   (A,C): 3
#   (B,C): 4
#   (A,B,C): 1

# Clean up
rm -f sample*.fq *.shared

# Generate fastq files
for m in 1 2; do
  ./random_fasta.py -reads 2 -len $bp > shared-AB_${m}.fq
  ./random_fasta.py -reads 3 -len $bp > shared-AC_${m}.fq
  ./random_fasta.py -reads 4 -len $bp > shared-BC_${m}.fq
  ./random_fasta.py -reads 1 -len $bp > shared-ABC_${m}.fq

  ./random_fasta.py -reads 4 -len $bp > sampleA_${m}.fq
  cat shared-AB_${m}.fq >> sampleA_${m}.fq
  cat shared-AC_${m}.fq >> sampleA_${m}.fq
  cat shared-ABC_${m}.fq >> sampleA_${m}.fq

  ./random_fasta.py -reads 3 -len $bp > sampleB_${m}.fq
  cat shared-AB_${m}.fq >> sampleB_${m}.fq
  cat shared-BC_${m}.fq >> sampleB_${m}.fq
  cat shared-ABC_${m}.fq >> sampleB_${m}.fq

  ./random_fasta.py -reads 2 -len $bp > sampleC_${m}.fq
  cat shared-AC_${m}.fq >> sampleC_${m}.fq
  cat shared-BC_${m}.fq >> sampleC_${m}.fq
  cat shared-ABC_${m}.fq >> sampleC_${m}.fq
done

# Generate the list of inputs, fqlist, that will be used for the -files argument
for i in A B C; do
  echo "sample${i} sample${i}_1.fq sample${i}_2.fq"
done | tr ' ' '\t' > fqlist

# First run only samples A and B
for i in A B; do
  fqsharedreads -files <(head -n 2 fqlist) -sample sample${i} > sample${i}.shared
done

# Next add sample C to the existing A and B runs.
fqsharedreads -continue sampleA.shared -files fqlist -sample sampleA > sampleA.continued.shared
fqsharedreads -continue sampleB.shared -files fqlist -sample sampleB > sampleB.continued.shared

# Also run sample C
fqsharedreads -files fqlist -sample sampleC > sampleC.shared
