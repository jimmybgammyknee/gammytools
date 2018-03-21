#!/bin/bash -l

# Jimmy Breen (jimmymbreen@gmail.com)
# 2015-02-05

# Call Variants using BCFtools from BAM file

if [ "$#" != "2" ]; then
        echo "Usage: $0 <reference> <sample_bam>"
        exit 0
fi 

samtools mpileup -R -ugf $1 $2 | bcftools call -vmO z -o "${basename $2 .bam}".vcf.gz
