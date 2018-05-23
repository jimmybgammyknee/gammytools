#!/bin/bash -l

# Usage: $0 [Data directory] [threads] [blacklist]

data=$1
threads=$2
blist=$3

base=$(pwd)
picard=$(which picard)
adjust=$(which tn5adjust.sh)

outdir=${base}/filter_clean

mkdir -p ${outdir}

for i in ${data}/_R1_exDupl.merged.bam; do
    samtools view -ubh -q 1 -F 2828 -@ ${threads} $i | \
            samtools sort -@ 10 -o ${outdir}/$(basename ${i} .merged.bam).filtered.bam -

    ${picard} MarkDuplicates REMOVE_DUPLICATES=TRUE \
        I=${outdir}/$(basename ${i} .merged.bam).filtered.bam \
        O=${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam \
        METRICS_FILE=${outdir}/$(basename ${i} .merged.bam).filtered.nodup.metrics.txt \
        VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE

    cat ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.metrics.txt | \
        cut -f 1 | grep -v chrM | \
            xargs samtools view -hb ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.bam > ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.bam

    bedtools subtract -A -abam ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.bam \
        -b ${blist} > ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.bam

    bash ${adjust} ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.bam \
        ${threads}

    samtools sort -@ 10 \
        -o ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.tn5Adj.sorted.bam \
            ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.tn5Adj.bam

    samtools index ${outdir}/$(basename ${i} .merged.bam).filtered.nodup.noMt.noblack.tn5Adj.sorted.bam
done
