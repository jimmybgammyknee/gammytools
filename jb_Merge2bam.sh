#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20141105
# Merge all PE samples in directory with bbmerge.sh and map merged reads to Reference Genome

if [ "$#" != "2" ]; then
        echo "Usage: jb_Merge2bam.sh [Sample directory] [BWA index]"
        exit 0
fi

#Load Modules
. /opt/shared/Modules/3.2.7/init/bash
module load picard bwa/0.7.9a 

IN=$1
REF=$2

cd $IN; for FQGZ in *_R1*.fastq.gz
 do
	/opt/local/bbmap/bbmerge.sh in1=$FQGZ in2=${FQGZ/R1/R2} out=${FQGZ%%.*}.merged.fastq.gz \
        	outu1=${FQGZ%%.*}.un1merged.fastq.gz outu2=${FQGZ%%.*}.un2merged.fastq.gz minlength=20 hist=${FQGZ%%.*}.merge.trim.hist
	bwa samse $REF <(bwa aln -t 32 $REF ${FQGZ%%.*}.merged.fastq.gz)  ${FQGZ%%.*}.merged.fastq.gz | \
		samtools view -Shb -q 30 /dev/stdin > HixHs37d5_${FQGZ%%.*}_BWA.bam
	java -jar /opt/shared/picard/1.71/SortSam.jar I=HixHs37d5_${FQGZ%%.*}_BWA.bam O=HixHs37d5_${FQGZ%%.*}_BWA.SORT.bam SO=coordinate
	java -jar /opt/shared/picard/1.71/MarkDuplicates.jar I=HixHs37d5_${FQGZ%%.*}_BWA.SORT.bam O=HixHs37d5_${FQGZ%%.*}_BWA.SORT.RMDUP.bam \
		AS=TRUE M=/dev/null REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT
done; cd ..
