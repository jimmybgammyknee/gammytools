#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20150205
# Call Variants using BCFtools from BAM file

if [ "$#" != "2" ]; then
        echo "Usage: jb_VCF_call_bcftools.sh [Reference] [BAM]"
        exit 0
fi 

ref=$1
bam=$2

bcf=/opt/local/samtools_1_1/bcftools-1.1/bcftools
sam=/opt/local/samtools_1_1/samtools-1.1/samtools

$sam mpileup -ugf $ref $bam | $bcf call -vmO z -o $bam.vcf.gz
