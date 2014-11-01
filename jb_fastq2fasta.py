#!/usr/bin/python

""" converts fastq to fasta on the fly """

import sys
from Bio import SeqIO
SeqIO.convert(sys.stdin, "fastq", sys.stdout, "fasta")
