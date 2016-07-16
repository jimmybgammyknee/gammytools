#!/bin/bash -l

## Jimmy Breen (jimmymbreen@gmail.com)

## Wrapper for the Salmon program
## $0 [R1 fastq] [R2 fastq] [index] [threads]

R1=$1
R2=$2
idx=$3
threads=$4

sn=$(basename $R1 .fastq.gz)

# Library stranded option is "IU"

salmon quant -i $idx \
	-l IU -p $threads \
	-1 <(pigz -cd $R1) -2 <(pigz -cd $R2) \
	-o ${sn}_${idx}_salmon
