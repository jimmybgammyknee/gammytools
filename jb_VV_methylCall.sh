#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com) -       20141118
# Calls methylation counts from fastafile and BAM (for the Grapevine Methylome paper)

if [ "$#" != "1" ]; then
        echo "Usage: jb_MethylCall.sh [BAM]"
        exit 0
fi

# I/O
BAM=$1
FastaRef=/rdsi/acad/Refs/BS_Genomes/Vvinifera_12X_chr/Vvinifera_12X_chr.fa
methylcall=/home/users/jbreen/bin/mpileup2methylation.py
python=/opt/shared/python/gcc4.4.4/2.7.2/bin/python
fregex=/home/users/jbreen/bin/fastaRegexFinder.py

cpg=/rdsi/acad/Refs/BS_Genomes/Vvinifera_12X_chr/Grapevine_BSseq.allcpg.bed.gz
chg=/rdsi/acad/Refs/BS_Genomes/Vvinifera_12X_chr/Grapevine_BSseq.allchg.bed.gz
chh=/rdsi/acad/Refs/BS_Genomes/Vvinifera_12X_chr/Grapevine_BSseq.allchh.bed.gz

# Loop through all bam files in the directory and call mC/C counts
samtools mpileup -l $cpg -f $FastaRef $BAM | \
	$python $methylcall -i - -f bismark | \
	awk '{print $1"\t"$2"\t"$3"\t""CG""\t"$4"\t"$5}' | pigz -c > ${BAM%%.*}.CG.txt.gz
samtools mpileup -l $chg -f $FastaRef $BAM | \
        $python $methylcall -i - -f bismark | \
        awk '{print $1"\t"$2"\t"$3"\t""CHG""\t"$4"\t"$5}' | pigz -c > ${BAM%%.*}.CHG.txt.gz
samtools mpileup -l $chh -f $FastaRef $BAM | \
        $python $methylcall -i - -f bismark | \
        awk '{print $1"\t"$2"\t"$3"\t""CHH""\t"$4"\t"$5}' | pigz -c > ${BAM%%.*}.CHH.txt.gz
