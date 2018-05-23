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
    macs2 callpeak --broad \            # MACS2 callpeaks function aggregating nearby peaks
        --format BAM \                      # let MACS2 use paired end information
        --nomodel --shift -37 --extsize 73 \  # Don't use default shift, instead use -37 and 73
        --bdg \                               # Create a bedgraph output
        --gsize 2.5e9 \                       # The mappable genome size
        -t $i \                               # Treatment file
        -n ${outdir}/$(basename ${i} .sorted.bam)        # Output prefix
done
