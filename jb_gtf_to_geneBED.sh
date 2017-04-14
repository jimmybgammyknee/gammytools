#!/bin/bash -l

# Converts a GTF file into a BED file with genes

# Jimmy Breen (jimmymbreen@gmail.com)

# jb_gtf_to_geneBED.sh [filein.gtf] [fileout.bed]

in=$1
out=$2

cat $in | \
    awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | \
    tr -d '";' > $out
