#!/bin/bash -l

# jimmymbreen@gmail.com	-	20141120
# Sort and Remove duplicates from Bismark Alignment bam

. /opt/shared/Modules/3.2.7/init/bash
module load samtools

samtools sort $1 ${basename $1 .bam}.sorted 
samtools rmdup -S ${basename $1 .bam}.sorted.bam ${basename $1 .bam}.sorted.noDup.bam
