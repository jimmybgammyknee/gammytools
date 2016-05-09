#!/bin/bash -l 

# Jimmy (jimmymbreen@gmail.com)	-	20141115
# Run AdapterRemoval on the R1 and R2 files

if [ "$#" != "3" ]; then
        echo "Usage: jb_AdapterRemoval_aDNA.sh [R1.fastq.gz] [R2.fastq.gz] [Adapters]"
        exit 0
fi

R1=$1
R2=$2
Adapters=$3

mkdir -p 2_Adapter

# Get adapters from the adapter file
ADAPT_5=$(head -n 1 $Adapters)
ADAPT_3=$(head -n 2 $Adapters | tail -n 1)

#Run
AdapterRemoval --file1 <(unpigz -c $R1) --file2 <(unpigz -c $R2) \
	--basename 2_Adapter/$(basename "${R1/_R1*/}") \
	--collapse --trimns --trimqualities \
	--minlength 25 --qualitybase 33 \
	--pcr1 $ADAPT_5 --pcr2 $ADAPT_3 --mm 3 \
	--output1 >(pigz --best > 2_Adapter/$(basename "${R1/_R1*/_1_truncated.fastq.gz}")) \
	--output2 >(pigz --best > 2_Adapter/$(basename "${R1/_R1*/_2_truncated.fastq.gz}")) \
	--discarded >(pigz --best > 2_Adapter/$(basename "${R1/_R1*/_discarded.fastq.gz}")) \
	--singleton >(pigz --best > 2_Adapter/$(basename "${R1/_R1*/_singletons.fastq.gz}")) \
	--outputcollapsed >(pigz --best > 2_Adapter/$(basename "${R1/_R1*/_collapsed.fastq.gz}")) \
	--outputcollapsedtruncated >(pigz --best > 2_Adapter/$(basename "${R1/_R1*/_collapsedTruncated.fastq.gz}"))
