#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20141002
# Process aligned BS-Seq data and produce methylation calls

if [ "$#" != "2" ]; then
        echo "Usage: jb_WGBS_bam2MethylCall.sh [Merged BAM file] [Reference Sequence file]"
        exit 0
fi

# Load modules
. /opt/shared/Modules/3.2.7/init/bash
module load bowtie/2-2.1.0
module load bismark
module load picard/1.71
module load python

# Inputs
BAM=$1
REF=$2

# Identify all CpG, CHG and CHH sites in the genome
fastaRegexFinder.py -f $REF -r CG --noreverse |pigz -c > ${REF%%.*}.allcpg.bed.gz 
fastaRegexFinder.py -f $REF -r C[ACT]G --noreverse |pigz -c > ${REF%%.*}.allchg.bed.gz
fastaRegexFinder.py -f $REF -r C[ACT]{2} --noreverse |pigz -c > ${REF%%.*}.allchh.bed.gz

# Run mpileup and create methylation calls 
samtools mpileup -l ${REF%%.*}.allcpg.bed.gz -f $REF $BAM | \
	 mpileup2methylation.py -i - -f bismark | pigz -c > ${BAM%%.*}.CpG.txt.gz
samtools mpileup -l ${REF%%.*}.allchg.bed.gz -f $REF $BAM | \
	 mpileup2methylation.py -i - -f bismark | pigz -c > ${BAM%%.*}.ChG.txt.gz
samtools mpileup -l ${REF%%.*}.allchh.bed.gz -f $REF $BAM | \
	 mpileup2methylation.py -i - -f bismark | pigz -c > ${BAM%%.*}.Chh.txt.gz
