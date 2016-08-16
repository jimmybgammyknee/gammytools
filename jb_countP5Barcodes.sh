#!/bin/bash -l

## Takes the R1 fastq file and prints the P5 barcodes present (if used). ##

fqin=$1

zcat $fqin | sed -n '2~4p' | cut -c1-5 | sort | uniq -c |sort -nr 
