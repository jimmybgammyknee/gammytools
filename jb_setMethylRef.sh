#!/bin/bash -l

# Create context files for a reference genome

if [ "$#" != "1" ]; then
        echo "Usage: jb_setRef.sh [Reference fasta]"
        exit 0
fi

#Load Modules
. /opt/shared/Modules/3.2.7/init/bash
module load bismark
module load bowtie/2-2.1.0
module load samtools
module load python/4.8.0/3.4.1

FastaRef=$1

# Progs
fregex=$(which fastaRegexFinder.py)
py34=/opt/shared/python/gcc4.8.0/3.4.1/bin/python3.4

$py34 $fregex -f $FastaRef -r CG --noreverse |pigz -c > ${FastaRef%%.*}.allcpg.bed.gz 
$py34 $fregex -f $FastaRef -r C[ACT]G --noreverse |pigz -c > ${FastaRef%%.*}.allchg.bed.gz
$py34 $fregex -f $FastaRef -r C[ACT]{2} --noreverse |pigz -c > ${FastaRef%%.*}.allchh.bed.gz
