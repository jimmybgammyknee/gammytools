#!/bin/bash -l

# Creates a distance matrix from multi-sample VCF file

vcf=$1
prefix=$2

plink=$(which plink)

$plink --vcf $vcf \
  --distance ibs \
  --allow-extra-chr \
  --out $prefix
