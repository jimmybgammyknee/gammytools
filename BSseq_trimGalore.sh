#!/bin/bash -l

# jimmymbreen@gmail.com	-	20141120
# Trim BSseq data with trim_galore (6bp 5' and adapters)

if [ "$#" != "2" ]; then
        echo "Usage: set_bsRef.sh [R1] [R2]"
        exit 0
fi

trimg=$(which trim_galore)

$trimg --paired --clip_R1 6 --clip_R2 6 -a AGATCGGAAGAGC -a2 CAAGCAGAAGACG --fastqc $1 $2
