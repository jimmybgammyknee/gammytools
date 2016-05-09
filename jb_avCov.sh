#!/bin/bash -l

bam=$1

samtools depth $bam  | \
	 awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
