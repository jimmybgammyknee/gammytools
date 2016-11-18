'''
Author: Jimmy Breen (jimmymbreen@gmail.com)
Date: 160917

Methylation pipeline for whole-methylome
- Tested on grapevine
'''

import time
import sys
import os
from parameters import *
from subprocess import PIPE, Popen, call

#Build bismark genome
def buildBismark(dir):
    """
    Build bismark formatted genome
    """
    # Check if dir exists
    os.makedirs(dir, exist_ok=True)

    # Use subprocess to run command
    subprocess('bismark_genome_preparation')

def bismarkRunPE(fq1, fq2, outdir):
    """
    Run bismark on paired-end data
    """

def bismarkRunSE(fq1, outdir):
    """
    Run bismark on single-end data
    """

def cleanUp(dirtybam, cleanbam):
    """
    Clean up bam file by sorting then removing duplicates
    """

    clean = subprocess.Popen("sambamba sort -p -t {0} {1} | sambamba markdup -p -r -t {0} /dev/stdin {2}".format(threads, dirtybam, cleanbam), shell=TRUE)
    return cleanbam
