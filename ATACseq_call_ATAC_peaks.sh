#!/bin/bash -l

# Call peaks on a sorted BAM for ATACseq analysis pipeline

# Jimmy Breen (jimmymbreen@gmail.com)
# 2018-05-23

# Input variable check
if [ $# -lt 1 ]; then
    echo "Usage: $0 [Data directory]"
    exit 1
fi

data=$1
base=$(pwd)

# Create output directory
outdir=${base}/macs2_peaks
mkdir -p ${outdir}

## Now we can run MACS2
for i in ${data}/*.sorted.bam; do
    macs2 callpeak --broad --format BAM \
        --nomodel --shift -37 --extsize 73 \
        --bdg --gsize 2.5e9 -t $i \
        -n ${outdir}/$(basename ${i} .sorted.bam)
done
