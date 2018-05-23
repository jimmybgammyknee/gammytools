#!/bin/bash -l

# Diff peaks ATACseq analysis pipeline with 3 reps

# Jimmy Breen (jimmymbreen@gmail.com)
# 2018-05-23

# Input variable check
if [ $# -lt 3 ]; then
    echo "Usage: $0 [Data directory]"
    exit 1
fi

S1=$1
S2=$2
S3=$3

# Create output directory
base=$(pwd)
outdir=${base}/macs2_diffPeaks_3rep
mkdir -p ${outdir}

macs2 bdgdiff -l 147 -g 73 \
    --t1 $(basename $S1)_treat_pileup.bdg \
    --t2 $(basename $S2)_treat_pileup.bdg \
    --c1 $(basename $S1)_control_lambda.bdg \
    --c2 $(basename $S2)_control_lambda.bdg \
    --o-prefix ${outdir}/2rep1

macs2 bdgdiff -l 147 -g 73 \
    --t1 $(basename $S1)_treat_pileup.bdg \
    --t2 $(basename $S3)_treat_pileup.bdg \
    --c1 $(basename $S1)_control_lambda.bdg \
    --c2 $(basename $S3)_control_lambda.bdg \
    --o-prefix ${outdir}/2rep2

macs2 bdgdiff -l 147 -g 73 \
    --t1 $(basename $S2)_treat_pileup.bdg \
    --t2 $(basename $S3)_treat_pileup.bdg \
    --c1 $(basename $S2)_control_lambda.bdg \
    --c2 $(basename $S3)_control_lambda.bdg \
    --o-prefix ${outdir}/2rep3

#wc -l *.bed > bdgdiff.txt        # count the number of peaks in bdgdiff output
#wc -l *.broadPeak > numPeaks.txt # count number of peaks in original macs2 calls
