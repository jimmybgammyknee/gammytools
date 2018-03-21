#!/bin/bash

# Jimmy Breen (jimmymbreen@gmail.com)

# ACAD Mapping protocol for the Ancient plants project
#       - Merges paired reads, maps with aDNA params and sorts/remove dups

if [ "$#" != "4" ]; then
        echo "Usage: ACAD_merge_map.sh FASTQ1 FASTQ2 ReferenceGenomeIndex num_of_threads"
        exit 0
fi

# Input variables and params
fq1=$1
fq2=$2
ref=$3
threads=$4

# programs
bbmerge=$(which bbmerge.sh)
bwa=$(which bwa)
picard=$(which picard)

${bbmerge} in1=$fq1 in2=$fq2 out=$(basename $fq1 .fastq.gz)_merged.fastq.gz

# Mapping merged read
${bwa} aln -t $threads -l 1024 -n 0.01 -o 2 \
        $ref $(basename $fq1 .fastq.gz)_merged.fastq.gz > $(basename $fq1 .fastq.gz)_BWA.sai

${bwa} samse $ref \
        $(basename $fq1 .fastq.gz)_BWA.sai $(basename $fq1 .fastq.gz)_merged.fastq.gz | \
                samtools view -q 30 -bSh -F4 -> $(basename $fq1 .fastq.gz)_merged_BWA.bam

# cleanup
rm "$(basename $fq1 .fastq.gz)"_BWA.sai

# Sorting and marking Dups
${picard} SortSam I=$(basename $fq1 .fastq.gz)_merged_BWA.bam \
        O=$(basename $fq1 .fastq.gz)_merged_BWA.SORT.bam \
        SO=coordinate

${picard} MarkDuplicates I=$(basename $fq1 .fastq.gz)_merged_BWA.SORT.bam \
        O=$(basename $fq1 .fastq.gz)_merged_BWA.SORT.RMDUP.bam \
        AS=TRUE M=/dev/null REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT

# cleanup
rm $(basename $fq1 .fastq.gz)_merged_BWA.bam
