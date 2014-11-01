#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com) 	-	20141101
# Take BAM file and produce bigWig annotation files

if [ "$#" != "2" ]; then
        echo "Usage: jb_bam2bigWig.sh [BAM file] [Reference]"
        exit 0
fi 

bamin=$1
ref=$2

cat $ref | python ~/bin/jb_fasta2lengthList.py > ${ref%%.*}.genome

bedtools genomecov -bg -split -strand + -ibam $bamin -g ${ref%%.*}.genome \
	| sort -k1,1 -k2,2n > ${bamin%%.*}.plus.bedGraph
bedtools genomecov -bg -split -strand - -ibam $bamin -g ${ref%%.*}.genome \
	| sort -k1,1 -k2,2n > ${bamin%%.*}.minus.bedGraph

bedGraphToBigWig ${bamin%%.*}.plus.bedGraph ${ref%%.*}.genome ${bamin%%.*}.plus.bigWig
bedGraphToBigWig ${bamin%%.*}.minus.bedGraph ${ref%%.*}.genome ${bamin%%.*}.minus.bigWig	
