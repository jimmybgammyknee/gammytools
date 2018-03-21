#!/bin/bash -l

# Jimmy Breen (jimmymbreen@gmail.com)
# 2018-03-21

# Methylation calling function using PileOMeth 
# 	- requires duplicate flagged BAM and reference file

# Input variables
file=$2
ref=$1

# Input variable checks
if [ $# -lt 2 ]
then 
	echo "Usage: $0 <ref> <bam>"
	exit 1
fi

# Index if none is present
if [ ! -e ${file}.bai ]	 
 then
	echo ">>> BAM file ${file} is being indexed"
	samtools index ${file}
fi

# run mbias function of PileOMeth to identify any strand bias 
PileOMeth mbias -q 30 --ignoreFlags --CHG --CHH \
	${ref} ${file} $(basename ${file} .bam) 2> flag.txt

# Read mbias flags into variable for extract function
opt=$(sed 's/^Suggested inclusion options://g' flag.txt)

PileOMeth extract -d 10 -q 30 ${opt} --CHG --CHH ${ref} ${file}

