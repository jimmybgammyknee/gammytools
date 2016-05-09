#!/usr/bin/env python

# Jimmy (jimmymbreen@gmail.com)	-	20141101
# Take Reference Genome and create .genome file for bedtools

import sys
from Bio import SeqIO

for record in SeqIO.parse(sys.stdin, "fasta"):
	print(record.id + "\t" + str(len(record.seq)))
