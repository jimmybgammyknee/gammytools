#!/bin/bash -l

sra_in=$1

module load ncbi/sratoolkit-2.2.2a

sam-dump $sra_in | samtools view -bS - | samtools sort - ${sra_in%%.*}.sorted && samtools index ${sra_in%%.*}.sorted.bam
