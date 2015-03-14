#!/bin/bash -l

# jimmymbreen@gmail.com	-	20141120
# Trim BSseq data with trim_galore (6bp 5' and adapters)

. /opt/shared/Modules/3.2.7/init/bash
module load fastQC
module load cutadapt

trimg='/opt/local/trim_galore_zip/trim_galore'
dir=$1

for FQGZ in ./$dir/*_R1*.fastq.gz
 do
	$trimg --paired --clip_R1 6 --clip_R2 6 -a AGATCGGAAGAGC -a2 CAAGCAGAAGACG --fastqc $FQGZ ${FQGZ/R1/R2}
done
