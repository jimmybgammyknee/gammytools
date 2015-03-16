#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com) -       20141118
# Calls methylation counts from fastafile and BAM (for the Grapevine Methylome paper)

if [ "$#" != "2" ]; then
        echo "Usage: jb_MethylCall.sh [BAM] [Reference fasta] [BS reference dir ]"
        exit 0
fi

. /opt/shared/Modules/3.2.7/init/bash
module load samtools

# I/O
bam=$1
faRef=$2
refDir=$3
mcall=/home/a1650598/bin/mpileup2methylation.py
py34=/opt/shared/python/gcc4.8.0/3.4.1/bin/python3.4
fregex=/home/a1650598/bin/fastaRegexFinder.py

cpg=$Refdir/${faRef%%.*}.allcpg.bed.gz.allcpg.bed.gz
chg=$Refdir/${faRef%%.*}.allchg.bed.gz.allchg.bed.gz
chh=$Refdir/${faRef%%.*}.allchh.bed.gz.allchh.bed.gz

# Loop through all bam files in the directory and call mC/C counts
samtools mpileup -l $cpg -f $faRef $bam | \
	$py34 $mcall -i - -f bismark | \
	awk '{print $1"\t"$2"\t"$3"\t""CG""\t"$4"\t"$5}' | pigz -c > ${bam%%.*}.CG.txt.gz
samtools mpileup -l $chg -f $faRef $bam | \
        $py34 $mcall -i - -f bismark | \
        awk '{print $1"\t"$2"\t"$3"\t""CHG""\t"$4"\t"$5}' | pigz -c > ${bam%%.*}.CHG.txt.gz
samtools mpileup -l $chh -f $faRef $bam | \
        $py34 $mcall -i - -f bismark | \
        awk '{print $1"\t"$2"\t"$3"\t""CHH""\t"$4"\t"$5}' | pigz -c > ${bam%%.*}.CHH.txt.gz
