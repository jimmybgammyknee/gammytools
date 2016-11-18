'''
Author: Jimmy Breen (jimmymbreen@gmail.com)
Date: 160917

Paramaters for all NGS processes
'''

#common
threads = 16

# Bismark
map_option = "--bowtie2"
score = "--score_min L,0,-0.6 -X 1000"
bam_out = "--bam"

#slurm

userMail = "jimmymbreen@gmail.com"
