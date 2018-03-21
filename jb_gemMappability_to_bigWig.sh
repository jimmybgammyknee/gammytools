#!/bin/bash -l

base=$(pwd)
REF=$1
THREADS=$2

# Load module
module load Biopython

# Create sizes file
~/bin/seq_length.py ${REF} > ${REF}.sizes

# Create index
gem-indexer -T ${THREADS} -c dna -i ${REF} -o $(basename ${REF}).gem

for kmer in $(seq 50 50 250); do

  # compute mappability data
  gem-mappability -T ${THREADS} -I $(basename ${REF}).gem \
    -l ${kmer} -o $(basename ${REF}).gem.${kmer}

  # convert results to wig and bigwig
  gem-2-wig -I $(basename ${REF}).gem \
    -i $(basename ${REF}).gem.${kmer}.mappability \
    -o $(basename ${REF}).gem.${kmer}
  wigToBigWig $(basename ${REF}).gem.${kmer}.wig \
    ${REF}.sizes \
    $(basename ${REF}).gem.${kmer}.bw

done
