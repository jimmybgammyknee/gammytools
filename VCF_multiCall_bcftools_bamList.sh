#!/bin/bash -l

# Jimmy Breen (jimmymbreen@gmail.com)
# 2015-02-05

# Call Variants in a group using BCFtools from all BAMs in list

if [ "$#" != "3" ]; then
        echo "Usage: $0 <reference> <bam_list> <vcf_prefix>"
        exit 0
fi 

ref=$1
bam=$2
prefix=$3

samtools mpileup -R -ugf ${ref} -b ${bam} | bcftools call -vmO z -o ${prefix}.vcf.gz
