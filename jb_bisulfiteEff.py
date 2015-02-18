#!/usr/bin/python env

# Estimate bisulfite conversion efficiency by counting the 
#	unconverted 5mC from mitochondrial or chloroplast sequences

# Bamfile must have the cp/mt sequence contained in the reference

import sys
import pysam as ps
import optparse

# function -> split bam into forward and reverse BAMs that 
#				match desired reference (CP or MT)

samfile = ps.AlignmentFile('./sdy10LV_bismark_bt2_pe.sorted.bam')
for pileupcolumn in samfile.pileup( 'chrCP' ):
	for read in pileupcolumn.pileups:
        	if pileupcolumn.nsegments == 10 and not read.indel and not read.is_del:
			print ("%s Base Coverge > 10 is actuially %s" % (pileupcolumn.pos, len(read.alignment))
		else:
			print ("Coverage < 10")
samfile.close()
