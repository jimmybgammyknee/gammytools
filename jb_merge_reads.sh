#!/bin/bash -l

# Need both java's to be stable on both machines
module load java
module load Java

in1=$1
in2=$2

bbmerge.sh in1=$1 \
           in2=$2 \
           out=$(basename $1 _R1.fastq.gz)_merged.fq.gz \
           outu1=$(basename $1 _R1.fastq.gz)_unmerged1.fq.gz \
           outu2=$(basename $1 _R1.fastq.gz)_unmerged1.fq.gz \
           outc=$(basename $1 _R1.fastq.gz)_kmer.txt \
           hist=$(basename $1 _R1.fastq.gz)_histogram.txt
