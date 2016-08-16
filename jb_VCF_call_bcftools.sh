#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20150205
# Call Variants using BCFtools from BAM file

if [ "$#" != "2" ]; then
        echo "Usage: jb_VCF_call_bcftools.sh [Reference] [BAM]"
        exit 0
fi 

. /opt/shared/Modules/3.2.7/init/bash
module load gnu/4.8.0
module load zlib/1.2.8-gnu_4.8.0
module load htslib/1.2.1
module load samtools/1.2 SamTools/1.2
module load bcftools/1.2


samtools mpileup -R -ugf $1 $2 | bcftools call -vmO z -o ${basename $1 .bam}_${basename $2}vcf.gz
