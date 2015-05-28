#!/bin/bash -l

# jimmymbreen@gmail.com
# Runs picard tools sort and markduplicates on single bam

bam=$1
name=$(basename $bam .bam)

java -jar /opt/shared/picard/1.71/SortSam.jar I=$bam O="$name".SORT.bam SO=coordinate
java -jar /opt/shared/picard/1.71/MarkDuplicates.jar I="$name".SORT.bam O="$name".SORT.RMDUP.bam AS=TRUE M=/dev/null REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT
