#!/bin/bash -l

# jimmymbreen@gmail.com	-	20141120
# Sort and Remove duplicates from Bismark Alignment bam

BAMIN=$1

samtools sort $BAMIN - | samtools rmdup -S - ${BAMIN%%.*}_sorted_noDup.bam
