#!/bin/bash -l

# Jimmy Breen (jimmymbreen@gmail.com)
# 2018-03-21

# Trim using trim galore, map with bwa-meth and 
# 	remove duplicates using Picard
#       - requires fastq files and reference file

ref=$1
fastq=$2
threads=$3

# Program variables
PICARD_HOME=/localscratch/Programs/bcbio/local/bin/picard
TRIMG=/localscratch/Programs/bcbio/anaconda/envs/wgbs/bin/trim_galore
BWA=/localscratch/Programs/bcbio/anaconda/envs/wgbs/bin/bwameth.py

# Input variable checks
if [ $# -lt 3 ]
then
        echo "Usage: $0 <ref> <fastq_directory> <threads>"
        exit 1
fi

# Define output directories
base=`pwd`
mkdir -p ${base}/1_trimGaloreOut
mkdir -p ${base}/2_bwaMethOut
mkdir 

# Index if none is present
if [ ! -e ${ref}.bwameth.c2t ]
 then
        echo ">>> Reference file ${ref} is being indexed"
	${BWA} index ${ref}
fi


# Loop processing over each samples
## NOTE: FILE EXTENSION ARE LIKELY TO CHANGE (e.g. _R1.fastq.gz or _1.fq.gz etc)

for FQGZ in ${fastq}/*_1*.fq.gz
 do

	echo ">>> Sample file ${FQGZ} trimmed using trim_galore and 6bp clip"

	# Trim galore
	${TRIMG} --paired --clip_R1 6 --clip_R2 6 -a AGATCGGAAGAGC -a2 CAAGCAGAAGACG \
        	--fastqc -o ${base}/1_trimGaloreOut/ ${FQGZ} ${FQGZ/_1/_2}

	echo ">>>  Sample file ${FQGZ} mapped using bwa-meth"

	# BWA-meth
	${BWA} --threads ${threads} --reference ${ref} \
		${base}/1_trimGaloreOut/$(basename ${FQGZ} _1.fq.gz)_1_val_1.fq.gz \
		${base}/1_trimGaloreOut/$(basename ${FQGZ} _1.fq.gz)_2_val_2.fq.gz | \
			samtools view -bhS - > ${base}/2_bwaMethOut/$(basename ${FQGZ} _1).bwameth_gatkBundle.bam

	# sambamba sort 
	sambamba sort ${base}/2_bwaMethOut/$(basename ${FQGZ} _1).bwameth_gatkBundle.bam
	
	echo ">>>  Sample file ${FQGZ} de-duplicate using Picard/MarkDuplicates"

	# Mark duplicates with picard and index output
	$PICARD_HOME MarkDuplicates INPUT=${base}/2_bwaMethOut/$(basename ${FQGZ} _1).bwameth_gatkBundle.sorted.bam \
        	OUTPUT=${base}/2_bwaMethOut/$(basename ${FQGZ} _1).bwameth_gatkBundle.sorted.markDups.bam \
       		METRICS_FILE=${base}/2_bwaMethOut/$(basename ${FQGZ} _1).bwameth_gatkBundle.sorted.markDups_metrics.txt \
        	REMOVE_DUPLICATES=false \
        	ASSUME_SORTED=true \
        	PROGRAM_RECORD_ID='null' \
        	VALIDATION_STRINGENCY=LENIENT
    	
	# Index bam
	samtools index ${base}/bwaMethOut/$(basename ${FQGZ} _1).bwameth_gatkBundle.sorted.markDups.bam

done
