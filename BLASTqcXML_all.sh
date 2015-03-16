#!/bin/bash

# Load Modules
. /opt/shared/Modules/3.2.7/init/bash
module load biopython
module load blast+

merged=$1
db=$2
threads=$3

if [ "$#" != "3" ]; then
        echo "Usage: BLASTQC2tbl.sh fastq.gz1 Database threads"
        exit 0
fi 

name=$(basename $merged .fastq.gz)
dbname=$(basename $db)

zcat $merged | sed -n '1~4s/^@/>/p;2~4p' | \
	/opt/shared/blast/blast-2.2.28+/bin/blastn -query - -db $db -num_threads $threads -evalue 1e-10 \
		-outfmt 5 -out "$name"_"$dbname".xml 
