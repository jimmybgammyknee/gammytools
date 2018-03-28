#!/bin/bash -l

# Jimmy Breeen (jimmymbreen@gmail.com)
# 2017-11-23

## Convert VCF to genotypes
##	- Must be gzipped
## $0 [filein]

file=$1

bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL ${file} | \
	cut -f1,2,4,5,6,12- | \
		pigz -c > $(basename $file .vcf.gz).genotypes.tsv.gz
