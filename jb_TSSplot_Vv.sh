#!/bin/bash -l

. /opt/shared/Modules/3.2.7/init/bash
module load gnu/4.8.0
module load zlib
module load bedtools

file=$1
tss=/localscratch/Refs/PlantGenomes/Vitis_vinifera/Vvinifera_145_Genoscope.12X.gene.gff3.tss.bed
sizes=/localscratch/Refs/PlantGenomes/Vitis_vinifera/Vvinifera_145_Genoscope_12X_Chromosomes.sizes

# extend tss
bedtools slop -b 5000 -i $tss -g $sizes > TSS_updown_5kb_slop.bed

# Prepare the tiled windows from TSS file
bedtools makewindows -b TSS_updown_5kb_slop.bed -w 50 -i winnum | sort -k1,1 -k2,2n > TSS_50bp_updown_5kb_slop.bed

# Make sure bed file is sorted
sort-bed <(zcat $file) > sorted.bed

# Get average methylation % over windows
# Get the average over ALL windows of all genes
bedtools map -a TSS_50bp_updown_5kb_slop.bed -b sorted.bed -c 5 -o mean -null 0 \
	| sort -t$'\t' -k4,4n \
	| bedtools groupby -i - -g 4 -c 5 -o mean > ${file}.TSS.counts.txt 

# Clean up
rm sorted.bed TSS_50bp_updown_5kb_slop.bed TSS_updown_5kb_slop.bed
