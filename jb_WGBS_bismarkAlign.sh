#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20141002
# Run Bismark on all Trimmed Sample in a directory in parallel
# (Whole genome BS-Seq data)

# Load modules
. /opt/shared/Modules/3.2.7/init/bash
module load bowtie/2-2.1.0
module load bismark
module load picard/1.71
module load python
module load parallel

if [ "$#" != "3" ]; then
        echo "Usage: jb_WGBS_align.sh [Data Directory] [Reference file directory] [num_of_threads]"
        exit 0
fi

if [ ! -d "$2" ]; then
    echo "Reference Genome \"$2\" doesnt exist. Now creating..."
    bismark_genome_preparation --bowtie2 $Refdir
    exit 1
else
    echo "Reference Genome is \"$2\" "
fi

# Inputs
data=$1
Refdir=$2

# Processes to run in parallel
threads=$3

# Run bismark alignment command in parallel
parallel -j $threads 'bismark -p 4 --bowtie2 --score_min L,0,-0.6 -X 1000 --bam \
	{1} -1 {2} -2 {3}' ::: $Refdir ::: $data/*_R1_001_val_1.fq.gz ::: $data/*_R2_001_val_2.fq.gz
wait 

parallel -j $threads 'samtools sort {} - | samtools rmdup -S - {}.sorted.nodup.bam' ::: *_bismark_bt2_pe.bam
