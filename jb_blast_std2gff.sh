#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com) -       20150328
# Produce gff file from a xml -> std tabular format

if [ "$#" != "1" ]; then
        echo "Usage: jb_blast_std2gff.sh [BLAST std tab]"
        exit 0
fi

stdin=$1

awk '{print $1"\t"".""\t""blast""\t"$7"\t"$8"\t"".""\t""+""\t"".""\t""Parent="$2}'  $stdin >  ${stdin%%.*}.gff
