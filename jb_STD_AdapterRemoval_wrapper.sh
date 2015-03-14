#!/bin/bash -l

# $0 [fastq R1] [fastq 2R] 
# AdapterRemoval with standard adapters

r1=$1
r2=$2

AdapterRemoval --file1 <(unpigz -c $r1) --file2 <(unpigz -c $r2) \
	--basename ${r1%%.*} --collapse --trimns --trimqualities \
	--minlength 25 --qualitybase 33 --mm 3 \
	--output1 >(pigz --best > ${r1%%.*}_1_truncated.fastq.gz) \
	--output2 >(pigz --best > ${r1%%.*}_2_truncated.fastq.gz) \
	--discarded >(pigz --best > ${r1%%.*}_discarded.fastq.gz) \
	--singleton >(pigz --best > ${r1%%.*}_singletons.fastq.gz) \
	--outputcollapsed >(pigz --best > ${r1%%.*}_collapsed.fastq.gz) \
	--outputcollapsedtruncated >(pigz --best > ${r1%%.*}_collapsedTruncated.fastq.gz)
