#!/bin/bash -l

# jimmymbreen@gmail.com	-	20141120
# Trim BSseq data with trim_galore (6bp 5' and adapters)

trimg='/opt/local/trim_galore_zip/trim_galore'

$trimg --paired --clip_R1 6 --clip_R2 6 -a AGATCGGAAGAGC -a2 CAAGCAGAAGACG --fastqc $R1 $R2
