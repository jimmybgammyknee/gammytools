#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com) -       20150205
# Call Variants using BCFtools from BAM file

if [ "$#" != "2" ]; then
        echo "Usage: jb_VCF_call_bcftools.sh [Reference] [bam_list]"
        exit 0
fi 

. /opt/shared/Modules/3.2.7/init/bash
module load gnu/4.8.0
module load zlib/1.2.8-gnu_4.8.0
module load htslib/1.2.1
module load samtools/1.2
module load bcftools/1.2

ref=$1
bam=$2

samtools mpileup -R -ugf $ref -b $bam | bcftools call -vmO z -o $bam.vcf.gz
