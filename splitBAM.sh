#!/bin/bash -l
# subsample BAM file for IDR on pseudoreplicates

SAMPLE=$1
THREADS=$2

NLINES=$( samtools view ${SAMPLE} | wc -l )
NLINES=$(( (NLINES + 1) / 2 ))

OUTNAME=${SAMPLE%.bam}_sample_

samtools view -H ${SAMPLE} > ${OUTNAME}00.sam
cp ${OUTNAME}00.sam ${OUTNAME}01.sam

samtools view ${SAMPLE} | shuf | split -d -l ${NLINES} - ${OUTNAME}
cat ${OUTNAME}00 >> ${OUTNAME}00.sam
cat ${OUTNAME}01 >> ${OUTNAME}01.sam

samtools view -hb -@ ${THREADS} ${OUTNAME}00.sam -o ${OUTNAME}00.bam
samtools view -hb -@ ${THREADS} ${OUTNAME}01.sam -o ${OUTNAME}01.bam

rm ${OUTNAME}00.sam ${OUTNAME}01.sam ${OUTNAME}00 ${OUTNAME}01
exit 0
