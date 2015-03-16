#!/bin/bash -l

# jimmymbreen@gmail.com	-	20141120
# Sort and Remove duplicates from Bismark Alignment bam

. /opt/shared/Modules/3.2.7/init/bash
module load samtools

BAMIN=$1

samtools sort $BAMIN ${BAMIN%%.*}_sorted 
samtools rmdup -S ${BAMIN%%.*}_sorted.bam ${BAMIN%%.*}_sorted_noDup.bam
