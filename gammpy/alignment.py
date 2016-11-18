'''
Alignment functions
'''

def bowtie2_alignPE(fq1, fq2, genome):
    cmd = ''

def mapped_and_sort(bam):
    samtools_view = ['samtools', 'view', '-bhS', '-']
    samtools_sort = ['samtools', 'sort', '-', prefix]
    samtools_index = ['samtools', 'index', bamfile]

    p1 = Popen(bowtie2_cmd, stdout = PIPE, stderr = bowtie2_logfh)
    p2 = Popen(samtools_view, stdin = p1.stdout, stdout = PIPE, stderr = bowtie2_logfh)
    p3 = Popen(samtools_sort, stdin = p2.stdout, stdout = PIPE, stderr = bowtie2_logfh)
    p1.stdout.close()
    p2.stdout.close()
