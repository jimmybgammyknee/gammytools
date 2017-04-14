#!/bin/bash -l

# usage: bam2cram.sh [Genome Reference] [BAM]

samb=$(which sambamba)
out="$(basename $2 .bam)".cram

$samb view -T $1 --format=cram -o $out $2
