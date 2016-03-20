#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20141002
# Merge WG BS-Seq replicates 

if [ "$#" != "3" ]; then
        echo "Usage: jb_WGBS_mergeReplicates.sh [rep1 BAM] [rep2 BAM] [rep3 BAM]"
        exit 0
fi

rep1=$1
rep2=$2
rep3=$3

samtools merge -nr ${rep1%%.*}.merged.bam $rep1 $rep2 $rep3