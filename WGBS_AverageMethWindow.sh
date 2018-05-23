#!/bin/bash -l

suff=$4
window=$3
dir=$2
ref=$1

# Input variable checks
if [ $# -lt 4 ]
then
        echo "Usage: $0 <ref> <directory containing bedGraph(gz)> <window_size> <file suffix>\n"
        echo "\n"
	echo "Example: $0 hg19.fasta directory 1000 _CpG.bedGraph.gz"
	exit 1
fi

mkdir -p ./"${window}"_windows_out
length=`which seq_length.py`

# Index if none is present
if [ ! -e ${ref}.genome ]
 then
        echo ">>> Genome file ${ref} is being created"
        ${length} ${ref} > ${ref}.genome
fi

# Make 1kb windows of the whole genome fasta
bedtools makewindows -g ${ref}.genome -w ${window} | awk -v w="${window}" '$3-$2 == w' > ./"${window}"_windows_out/$(basename ${ref})."${window}".windows

# Loops through all WGBS methylation bedGraph files and produce 1kb tiles of average methylation 
for bg in ${dir}/*"${suff}"
 do
	sortBed -i ${bg} | bedtools map -a ./"${window}"_windows_out/$(basename ${ref})."${window}".windows -b stdin -c 4 -o mean > ./"${window}"_windows_out/${bg}.amc."${window}".windows
done	

# Merge all bedGraph 1kb window tiles and remove empty intervals
bedtools unionbedg -header -i ./"${window}"_windows_out/*.amc."${window}".windows > ./"${window}"_windows_out/union.amc."${window}".windows
