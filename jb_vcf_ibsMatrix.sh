#!/bin/bash -l

# Creates a distance matrix from multi-sample VCF file

vcf=$1
prefix=$2

plink=$(which plink1.90)

$plink --vcf $vcf \
  --distance ibs \
  --allow-extra-chr \
  --out $prefix
