#!/bin/bash -l

# BAM filtering script for ATACseq analysis pipeline

# Jimmy Breen (jimmymbreen@gmail.com)
# 2018-05-23

# Input variable check
if [ $# -lt 3 ]; then
    echo "Usage: $0 [Data directory] [threads] [blacklist]"
    exit 1
fi

data=$1
threads=$2
blist=$3

base=$(pwd)
picard=$(which picard)
adjust=$(which tn5adjust.sh)

# Create output directory
outdir=${base}/filter_clean
mkdir -p ${outdir}

# Loop through each bam file in data directory
for i in ${data}/*.bam; do

    # Discard unmapped, non-primary, qc failed, and supp alignments
    samtools view -ubh -q 1 -F 2828 -@ ${threads} $i | \
            samtools sort -@ 10 -o ${outdir}/$(basename ${i} .merged.bam).filtered.bam -

    # Remove duplicates
    ${picard} MarkDuplicates REMOVE_DUPLICATES=TRUE \
        I=${outdir}/$(basename ${i} .merged.bam).filtered.bam \
        O=${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam \
        METRICS_FILE=${outdir}/$(basename ${i} .merged.bam).filtered.nodup.metrics.txt \
        VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE

    # Get index stats for mtDNA removal
    samtools index ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam
    samtools idxstats ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam > ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam.idxStats.txt

    cat ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam.idxStats.txt | \
        cut -f 1 | grep -v chrM | \
            xargs samtools view -hb ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam > ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.bam

    # Get rid of blacklisted sites from hg19
    bedtools subtract -A -abam ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.bam \
        -b ${blist} > ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.bam

    # Because the Tn5 transposase cuts in a staggered fashion, the resulting library of reads have + strand shifted downstream 4bp while the - strand is shifted upstream 5bp. Reversing this improves the accuracy of the open chromatin prediction.
    bash ${adjust} ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.bam \
        ${threads}

    # Sort and reindex
    samtools sort -@ 10 \
        -o ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.tn5Adj.sorted.bam \
            ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.tn5Adj.bam
    samtools index ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.tn5Adj.sorted.bam

done
