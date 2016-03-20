#!/bin/bash -l 

# Load required modules for R
. /opt/shared/Modules/3.2.7/init/bash
module load gnu/4.8.0 java/java-jdk-1.7.051
module load R/3.1.1
module load python/4.8.0/3.4.1
module load parallel

# IO
htseq=/opt/local/HTSeq-0.6.1/scripts/htseq-count
bamdir=$1
gff3=$2

# Run htseq-count on all bams
for i in $bamdir/*.bam
 do
	samtools view -h $i > ${i}.sam 
	python $htseq -i ID -t gene ${i}.sam $gff3 > ${i}.counts.txt  
    	rm $i.sam
done
