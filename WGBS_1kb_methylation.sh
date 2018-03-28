#!/bin/bash -l

suff=$3
dir=$2
ref=$1

# Input variable checks
if [ $# -lt 3 ]
then
        echo "Usage: $0 <ref> <directory containing bedGraph(gz)> <file suffix>\n"
        echo "\n"
	echo "Example: $0 hg19.fasta directory _CpG.bedGraph.gz"
	exit 1
fi

length=`which seq_length.py`

# Index if none is present
if [ ! -e ${ref}.genome ]
 then
        echo ">>> Genome file ${ref} is being created"
        ${length} ${ref} > ${ref}.genome
fi

# Make 1kb windows of the whole genome fasta
bedtools makewindows -g ${ref}.genome -w 1000 > ./$(basename ${ref}).1kb.windows

# Loops through all WGBS methylation bedGraph files and produce 1kb tiles of average methylation 
for bg in ${dir}/*"${suff}"
 do
	bedtools map -a ./$(basename ${ref}).1kb.windows -b ${bg} -c 4 -o mean > ${bg}.amc.1kb.windows
done	

# Merge all bedGraph 1kb window tiles and remove empty intervals
bedtools unionbedg -header -i *.amc.1kb.windows > union.amc.1kb.windows
