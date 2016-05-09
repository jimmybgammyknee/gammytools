#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com) -       20141118
# Calls methylation counts from fastafile and BAM (for the Grapevine Methylome paper)

if [ "$#" != "2" ]; then
        echo "Usage: jb_MethylCall.sh [REF] [BAM]"
        exit 0
fi

methylcall=$(which mpileup2methylation.py)
fregex=$(which fastaRegexFinder.py)
py34=$(which python3.4)

$py34 $fregex -f $1 -r CG --noreverse |pigz -c > ${basename $1}.allcpg.bed.gz
$py34 $fregex -f $1 -r C[ACT]G --noreverse |pigz -c > ${basename $1}.allchg.bed.gz
$py34 $fregex -f $1 -r C[ACT]{2} --noreverse |pigz -c > ${basename $1}.allchh.bed.gz

# Loop through all bam files in the directory and call mC/C counts
samtools mpileup -l ${basename $1}.allcpg.bed.gz -f $1 $2 | \
	python $methylcall -i - -f bismark | \
	awk '{print $2"\t"$2"\t"$3"\t""CG""\t"$4"\t"$5}' | pigz -c > ${basename $2}.CG.txt.gz
samtools mpileup -l ${basename $1}.allchg.bed.gz -f $1 $2 | \
        python $methylcall -i - -f bismark | \
        awk '{print $2"\t"$2"\t"$3"\t""CHG""\t"$4"\t"$5}' | pigz -c > ${basename $2}.CHG.txt.gz
samtools mpileup -l ${basename $1}.allchh.bed.gz -f $1 $2 | \
        python $methylcall -i - -f bismark | \
        awk '{print $2"\t"$2"\t"$3"\t""CHH""\t"$4"\t"$5}' | pigz -c > ${basename $2}.CHH.txt.gz
