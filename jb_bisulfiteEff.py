#!/usr/bin/python env

# Estimate bisulfite conversion efficiency by counting the 
#	unconverted 5mC from mitochondrial or chloroplast sequences

# Bamfile must have the cp/mt sequence contained in the reference

import sys
import pysam as ps
import optparse

import pyximport 
pyximport.install()
import _pysam_flagstat

# function -> split bam into forward and reverse BAMs that 
#				match desired reference (CP or MT)

samfile = ps.AlignmentFile('./sdy10LV_bismark_bt2_pe.sorted.bam')
cdef AlignedSegment read
for read in samfile:
	for pileupcolumn in samfile.pileup( 'chrCP' ):
		for read in pileupcolumn.pileups:
        		if pileupcolumn.nsegments == 10 and not read.indel and not read.is_del:
				if read.is_reverse == False:
					Cbase = (read.query_sequence.count('T') + read.query_sequence.count(','))
					print Cbase
					#Cbfreq = (Cbase * 1.0 / len(read.query_sequence)
					#print  pileupcolumn.pos#, Cbfreq
				elif read.is_reverse == True:
					Tbase = (read.query_sequence.count('a') + read.query_sequence.count('.'))
					print Tbase
					#Tbfreq = (Cbase * 1.0 / len(read.query_sequence)
					print pileupcolumn.pos#, Tbfreq
samfile.close()

# Grab only ref seq using ID
def reffetch(refid):
	"""
	reffetch - returns all chloroplast/mtDNA pileups using a reference ID
		e.g. 'chrCP'
	"""

# Feed in agruments
def GetArgs():
	"""
	GetArgs - read the command line
	returns - a bam_file and output file, genome file names
	
	typing python bisulfiteEff.py will show the options
	"""
	
	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
		def error(self, message):
			sys.stderr.write('error: %s\n' % message)
			self.print_help()
			sys.exit(2)

        parser = Parser(description='Calculate the bisulfite conversion efficiency of bam file')
        parser.add_argument('-b', '--bam_file',
			type=str,
			required=True,
			help='Bam file (required).')
        parser.add_argument('-o', '--output_file',
			type=str,
			help='output csv file.')
        parser.add_argument('-g', '--genome_file',
			required=True,
			type=str,
			help='File containing entire reference genome.')
	parser.add_argument('-i', '--reference_id',
			required=True,
			type=str,
			help='The Reference identifier for the sequence you want to use for calculations')

	return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.bam_file, args.output_file, args.genome_file, args.ref_id
