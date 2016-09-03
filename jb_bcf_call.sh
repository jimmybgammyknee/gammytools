#!/bin/bash -l

# Variant calling using bcftools - Needs a reference and bam list and prefix

if [[ $# -eq 0 ]] ; then
    echo 'Usage: jb_bcf_call.sh [bam list] [reference fasta] [outfile prefix]'
    exit 1
fi

samtools=/localscratch/Programs/samtools-1.2/bin/samtools
bcftools=/localscratch/Programs/bcftools-1.2/bin/bcftools
bam=$1
ref=$2
prefix=$3

$samtools mpileup -m 2 -p -F 0.1 -ugf $ref -b $bam | \
	$bcftools call -vmO v | \
	$bcftools norm -f $ref | \
	$bcftools filter -i "QUAL > 30 && DP > 10" | \
	bgzip > $prefix.samtools.norm.flt.vcf.gz
