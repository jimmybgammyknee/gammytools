#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20141002
# Run Bismark on Whole genome BS-Seq data

if [ "$#" != "2" ]; then
        echo "Usage: jb_WGBS_align.sh [Data Directory] [Reference file directory]"
        exit 0
fi

# Load modules
. /opt/shared/Modules/3.2.7/init/bash
module load bowtie/2-2.1.0
module load bismark
module load picard/1.71
module load python

# Inputs
data=$1
Refdir=$2
#threads=$3

# Prep genome
bismark_genome_preparation --bowtie2 $Refdir

# Loop over all R1/R2 files and align to BS Genome
for FQGZ in $data/*_R1*.fastq.gz
 do
        bismark -p 4 --bowtie2 --score_min L,0,-0.6 -X 1000 --bam $Refdir -1 $FQGZ -2 ${FQGZ/R1/R2}
        samtools sort ${FQGZ}_bismark_bt2_pe.bam - | samtools rmdup -S - ${FQGZ}_bismark_bt2_pe.sorted.nodup.bam
done