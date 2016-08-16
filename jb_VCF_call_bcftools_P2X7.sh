#!/bin/bash -l

. /opt/shared/Modules/3.2.7/init/bash
module load samtools/1.2 SamTools/1.2
module load zlib gnu/4.8.0
module load bcftools/1.2

ref=$1
bamin=$2

samtools index $bamin
samtools mpileup -ugf $ref -r 12:121570622-121623876 $bamin | bcftools call -vmO z -o ${bamin%%.*}.p2x7.vcf.gz
