#!/bin/bash -l

. /opt/shared/Modules/3.2.7/init/bash
module load gnu/4.8.0
module load zlib
module load bedtools

# Use annotation in full GFF3 annotation and make a TSS plot
# $0 [bedGraph.gz] [Annotation Type] [full gff3] [sizes file] [flanking]

file=$1
full=$3
sizes=$4
flank=$5

tag="$2"
zgrep "$tag" $full | awk '$7=="+" {print $1"\t"$4"\t"$4+1"\t"$7"\t"$9}' > "$tag"_start.bed
zgrep "$tag" $full | awk '$7=="-" {print $1"\t"$5-1"\t"$5"\t"$7"\t"$9}' >> "$tag"_start.bed

# Remove description section - not that it makes much of a difference
#sed -e 's/^/chr/g' -e 's/assembly_name=IGGP_12x;description=//g' "$tag"_start.bed
sed -e 's/^/chr/g' "$tag"_start.bed | \
	awk -v fl="$flank" '$2 > fl {print $0}' | \
	sort -k1,1 -k2,2n > "$tag"_start_sorted.bed 

echo "TSS file done"

# extend tss to 2kb
bedtools slop -b $5 -i "$tag"_start_sorted.bed -g $sizes > TSS_updown_window.bed

# Prepare the tiled windows from TSS file
bedtools makewindows -b TSS_updown_window.bed -w 50 -i winnum | \
	sort -k1,1 -k2,2n -k3,3n | sed '/^81/d' > TSS_50bp_updown_window.bed

echo "TSS window file done"

# Make sure bed file is sorted
sortBed -i <(zcat $file) > sorted.bed

echo "Sorted bed file"

# Get average methylation % over windows
# Get the average over ALL windows of all genes
bedtools map -a TSS_50bp_updown_window.bed -b sorted.bed -c 5 -o mean -null 0 \
        | sort -t$'\t' -k4,4n \
        | bedtools groupby -i - -g 4 -c 5 -o mean > ${file}.TSS.counts.txt

# Clean up
rm sorted.bed TSS_50bp_updown_window.bed TSS_updown_window.bed
