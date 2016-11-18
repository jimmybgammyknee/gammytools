'''
Author: Jimmy Breen (jimmymbreen@gmail.com)
Date:160918

Utilities for pipelines etc

Some really good functions borrowed from
https://github.com/afrendeiro
'''


def slurmHeader(jobName, output, queue="shortq", ntasks=1, time="10:00:00",
                cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
    cmd = """    #!/bin/bash
    #SBATCH --partition={0}
    #SBATCH --ntasks={1}
    #SBATCH --time={2}
    #SBATCH --cpus-per-task={3}
    #SBATCH --mem-per-cpu={4}
    #SBATCH --nodes={5}
    #SBATCH --job-name={6}
    #SBATCH --output={7}
    #SBATCH --mail-type=end
    #SBATCH --mail-user={8}
    # Start running the job
    hostname
    date
    """.format(queue, ntasks, time, cpusPerTask, memPerCpu,
               nodes, jobName, output, userMail)

    return cmd


def slurmFooter():
    cmd = "    date"

    return cmd


def slurmSubmitJob(jobFile):
    import os
    cmd = "sbatch %s" % jobFile

    return os.system(cmd)


def removeFile(fileName):
    cmd = "rm {0}".format(fileName)

    return cmd


def moveFile(old, new):
    cmd = "mv {0} {1}".format(old, new)

    return cmd


def makeDir(directory):
    cmd = "mkdir -p {0}".format(directory)

    return cmd


def mergeBams(inputBams, outputBam):
    cmd = "java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/MergeSamFiles.jar"
    cmd += " USE_THREADING=TRUE"
    cmd += " " + (" ".join(["INPUT=%s"] * len(inputBams))) % tuple(inputBams)
    cmd += " OUTPUT={0}".format(outputBam)

    return cmd


def fastqc(inputBam, outputDir, sampleName):
    import os

    initial = os.path.splitext(os.path.basename(inputBam))[0]

    cmd1 = "fastqc --noextract --outdir {0} {1}".format(outputDir, inputBam)

    cmd2 = "mv {0}_fastqc.html {1}_fastqc.html".format(os.path.join(outputDir, initial), os.path.join(outputDir, sampleName))

    cmd3 = "mv {0}_fastqc.zip {1}_fastqc.zip".format(os.path.join(outputDir, initial), os.path.join(outputDir, sampleName))

    return [cmd1, cmd2, cmd3]


def bam2fastq(inputBam, outputFastq, outputFastq2=None, unpairedFastq=None):
    cmd = "java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/SamToFastq.jar"
    cmd += " INPUT={0}".format(inputBam)
    cmd += " FASTQ={0}".format(outputFastq)
    if outputFastq2 is not None and unpairedFastq is not None:
        cmd += " SECOND_END_FASTQ={0}".format(outputFastq2)
        cmd += " UNPAIRED_FASTQ={0}".format(unpairedFastq)

    return cmd


def trimmomatic(inputFastq1, outputFastq1, cpus, adapters, log,
                inputFastq2=None, outputFastq1unpaired=None,
                outputFastq2=None, outputFastq2unpaired=None):

    PE = False if inputFastq2 is None else True
    pe = "PE" if PE else "SE"

    cmd = "java -Xmx4g -jar `which trimmomatic-0.32.jar`"
    cmd += " {0} -threads {1} -trimlog {2} {3}".format(pe, cpus, log, inputFastq1)
    if PE:
        cmd += " {0}".format(inputFastq2)
    cmd += " {0}".format(outputFastq1)
    if PE:
        cmd += " {0} {1} {2}".format(outputFastq1unpaired, outputFastq2, outputFastq2unpaired)
    cmd += " ILLUMINACLIP:{0}:1:40:15:8:true".format(adapters)
    cmd += " LEADING:3 TRAILING:3"
    cmd += " SLIDINGWINDOW:4:10"
    cmd += " MINLEN:36"

    return cmd


def bowtie2Map(inputFastq1, outputBam, log, metrics, genomeIndex, maxInsert, cpus, inputFastq2=None):
    import re

    outputBam = re.sub("\.bam$", "", outputBam)
    # Admits 2000bp-long fragments (--maxins option)
    cmd = "bowtie2 --very-sensitive -p {0}".format(cpus)
    cmd += " -x {0}".format(genomeIndex)
    cmd += " --met-file {0}".format(metrics)
    if inputFastq2 is None:
        cmd += " {0} ".format(inputFastq1)
    else:
        cmd += " --maxins {0}".format(maxInsert)
        cmd += " -1 {0}".format(inputFastq1)
        cmd += " -2 {0}".format(inputFastq2)
    cmd += " 2> {0} | samtools view -S -b - | samtools sort - {1}".format(log, outputBam)

    return cmd


def topHatMap(inputFastq, outDir, genome, transcriptome, cpus):
    # TODO:
    # Allow paired input
    cmd = "tophat --GTF {0} --b2-L 15 --library-type fr-unstranded --mate-inner-dist 120".format(transcriptome)
    cmd += " --max-multihits 100 --no-coverage-search"
    cmd += " --num-threads {0} --output-dir {1} {2} {3}".format(cpus, outDir, genome, inputFastq)

    return cmd

def picardMarkDuplicates(inputBam, outputBam, metricsFile, tempDir="."):
    import re

    transientFile = re.sub("\.bam$", "", outputBam) + ".dups.nosort.bam"
    outputBam = re.sub("\.bam$", "", outputBam)

    cmd1 = "java -Xmx4g -jar `which MarkDuplicates.jar`"
    cmd1 += " INPUT={0}".format(inputBam)
    cmd1 += " OUTPUT={0}".format(transientFile)
    cmd1 += " METRICS_FILE={0}".format(metricsFile)
    cmd1 += " VALIDATION_STRINGENCY=LENIENT"
    cmd1 += " TMP_DIR={0}".format(tempDir)

    # Sort bam file with marked duplicates
    cmd2 = "samtools sort {0} {1}".format(transientFile, outputBam)

    # Remove transient file
    cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(transientFile)

    return [cmd1, cmd2, cmd3]


def removeDuplicates(inputBam, outputBam, cpus=16):
    cmd = "sambamba markdup -t {0} -r {1} {2}".format(cpus, inputBam, outputBam)

    return cmd


def filterReads(inputBam, outputBam, metricsFile, paired=False, cpus=16, Q=30):
    """
    Remove duplicates, filter for >Q, remove multiple mapping reads.
    For paired-end reads, keep only proper pairs.
    """
    import re

    nodups = re.sub("\.bam$", "", outputBam) + ".nodups.nofilter.bam"

    cmd1 = "sambamba markdup -t {0} -r --compression-level=0 {1} {2} 2> {3}".format(cpus, inputBam, nodups, metricsFile)

    cmd2 = ' sambamba view -t {0} -f bam --valid'.format(cpus)
    if paired:
        cmd2 += ' -F "not (unmapped or mate_is_unmapped) and proper_pair'
    else:
        cmd2 += ' -F "not unmapped'
    cmd2 += ' and not (secondary_alignment or supplementary) and mapping_quality >= {0}"'.format(Q)
    cmd2 += ' {0} |'.format(nodups)
    cmd2 += " sambamba sort -t {0} /dev/stdin -o {1}".format(cpus, outputBam)

    cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(nodups)
    cmd4 = "if [[ -s {0} ]]; then rm {0}; fi".format(nodups + ".bai")

    return [cmd1, cmd2, cmd3, cmd4]


def sortIndexBam(inputBam, outputBam):
    import re

    tmpBam = re.sub("\.bam", ".sorted", inputBam)

    cmd1 = "samtools sort {0} {1}".format(inputBam, tmpBam)

    cmd2 = "mv {0}.bam {1}".format(tmpBam, outputBam)

    cmd3 = "samtools index {0}".format(outputBam)

    return [cmd1, cmd2, cmd3]


def indexBam(inputBam):
    cmd = "samtools index {0}".format(inputBam)

    return cmd


def getFragmentSizes(bam, maxInsert=1500):
    try:
        import pysam
        import numpy as np
    except:
        return

    fragSizes = list()

    bam = pysam.Samfile(bam, 'rb')

    for read in bam:
        if bam.getrname(read.tid) != "chrM" and read.tlen < maxInsert:
            fragSizes.append(read.tlen)
    bam.close()

    return np.array(fragSizes)


def plotInsertSizesFit(bam, plot, outputCSV, maxInsert=1500, smallestInsert=30):
    """
    Heavy inspiration from here:
    https://github.com/dbrg77/ATAC/blob/master/ATAC_seq_read_length_curve_fitting.ipynb
    """
    try:
        import pysam
        import numpy as np
        import matplotlib.mlab as mlab
        from scipy.optimize import curve_fit
        from scipy.integrate import simps
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except:
        print("Necessary Python modules couldn't be loaded.")
        return

    try:
        import seaborn as sns
        sns.set_style("whitegrid")
    except:
        pass

    def getFragmentSizes(bam, maxInsert=1500):
        fragSizes = list()

        bam = pysam.Samfile(bam, 'rb')

        for i, read in enumerate(bam):
            if read.tlen < maxInsert:
                fragSizes.append(read.tlen)
        bam.close()

        return np.array(fragSizes)

    def mixtureFunction(x, *p):
        """
        Mixture function to model four gaussian (nucleosomal) and one exponential (nucleosome-free) distributions.
        """
        m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = p
        nfr = expo(x, 2.9e-02, 2.8e-02)
        nfr[:smallestInsert] = 0

        return (mlab.normpdf(x, m1, s1) * w1 +
                mlab.normpdf(x, m2, s2) * w2 +
                mlab.normpdf(x, m3, s3) * w3 +
                mlab.normpdf(x, m4, s4) * w4 +
                nfr)

    def expo(x, q, r):
        """
        Exponential function.
        """
        return q * np.exp(-r * x)

    # get fragment sizes
    fragSizes = getFragmentSizes(bam)

    # bin
    numBins = np.linspace(0, maxInsert, maxInsert + 1)
    y, scatter_x = np.histogram(fragSizes, numBins, density=1)
    # get the mid-point of each bin
    x = (scatter_x[:-1] + scatter_x[1:]) / 2

    # Parameters are empirical, need to check
    paramGuess = [
        200, 50, 0.7,  # gaussians
        400, 50, 0.15,
        600, 50, 0.1,
        800, 55, 0.045,
        2.9e-02, 2.8e-02  # exponential
    ]

    try:
        popt3, pcov3 = curve_fit(mixtureFunction, x[smallestInsert:], y[smallestInsert:], p0=paramGuess, maxfev=100000)
    except:
        print("Nucleosomal fit could not be found.")
        return

    m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = popt3

    # Plot
    plt.figure(figsize=(12, 12))

    # Plot distribution
    plt.hist(fragSizes, numBins, histtype="step", ec="k", normed=1, alpha=0.5)

    # Plot nucleosomal fits
    plt.plot(x, mlab.normpdf(x, m1, s1) * w1, 'r-', lw=1.5, label="1st nucleosome")
    plt.plot(x, mlab.normpdf(x, m2, s2) * w2, 'g-', lw=1.5, label="2nd nucleosome")
    plt.plot(x, mlab.normpdf(x, m3, s3) * w3, 'b-', lw=1.5, label="3rd nucleosome")
    plt.plot(x, mlab.normpdf(x, m4, s4) * w4, 'c-', lw=1.5, label="4th nucleosome")

    # Plot nucleosome-free fit
    nfr = expo(x, 2.9e-02, 2.8e-02)
    nfr[:smallestInsert] = 0
    plt.plot(x, nfr, 'k-', lw=1.5, label="nucleosome-free")

    # Plot sum of fits
    ys = mixtureFunction(x, *popt3)
    plt.plot(x, ys, 'k--', lw=3.5, label="fit sum")

    plt.legend()
    plt.xlabel("Fragment size (bp)")
    plt.ylabel("Density")
    plt.savefig(plot, bbox_inches="tight")

    # Integrate curves and get areas under curve
    areas = [
        ["fraction", "area under curve", "max density"],
        ["Nucleosome-free fragments", simps(nfr), max(nfr)],
        ["1st nucleosome", simps(mlab.normpdf(x, m1, s1) * w1), max(mlab.normpdf(x, m1, s1) * w1)],
        ["2nd nucleosome", simps(mlab.normpdf(x, m2, s2) * w1), max(mlab.normpdf(x, m2, s2) * w2)],
        ["3rd nucleosome", simps(mlab.normpdf(x, m3, s3) * w1), max(mlab.normpdf(x, m3, s3) * w3)],
        ["4th nucleosome", simps(mlab.normpdf(x, m4, s4) * w1), max(mlab.normpdf(x, m4, s4) * w4)]
    ]

    try:
        import csv

        with open(outputCSV, "w") as f:
            writer = csv.writer(f)
            writer.writerows(areas)
    except:
        pass


def bamToBigWig(inputBam, outputBigWig, genomeSizes, genome, tagmented=False, normalize=False):
    import os
    import re
    import subprocess

    # TODO:
    # addjust fragment length dependent on read size and real fragment size
    # (right now it asssumes 50bp reads with 180bp fragments)

    cmds = list()

    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))

    cmd1 = "bedtools bamtobed -i {0} |".format(inputBam)
    if not tagmented:
        cmd1 += " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
        cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
    cmd1 += " genomeCoverageBed {0}-bg -g {1} -i stdin > {2}.cov".format(
            "-5 " if tagmented else "",
            genomeSizes,
            transientFile
    )
    cmds.append(cmd1)

    if normalize:
        cmds.append("""awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov""".format(transientFile))

    cmds.append("bedGraphToBigWig {0}{1}.cov {2} {3}".format(transientFile, ".normalized" if normalize else "", genomeSizes, outputBigWig))

    # remove tmp files
    cmds.append("if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transientFile))
    if normalize:
        cmds.append("if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(transientFile))

    cmds.append("chmod 755 {0}".format(outputBigWig))

    return cmds


def addTrackToHub(sampleName, trackURL, trackHub, colour, fivePrime=""):
    cmd1 = """echo "track type=bigWig name='{0} {1}' description='{0} {1}'""".format(sampleName, fivePrime)
    cmd1 += """ height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}" >> {2}""".format(trackURL, colour, trackHub)

    cmd2 = "chmod 755 {0}".format(trackHub)

    return [cmd1, cmd2]


def linkToTrackHub(trackHubURL, fileName, genome):
    import textwrap

    db = "org" if genome == "hg19" else "db"  # different database call for human
    genome = "human" if genome == "hg19" else genome  # change hg19 to human

    html = """
    <html>
        <head>
            <meta http-equiv="refresh" content="0; url=http://genome.ucsc.edu/cgi-bin/hgTracks?"""
    html += """{db}={genome}&hgt.customText={trackHubURL}" />
        </head>
    </html>
    """.format(trackHubURL=trackHubURL, genome=genome, db=db)

    with open(fileName, 'w') as handle:
        handle.write(textwrap.dedent(html))


def htSeqCount(bam, gtf, output):

    #Changed for bams
    cmd = "htseq-count -f bam -t gene -i gene_id {0} {1} > {2}".format(bam, gtf, output)

    return [cmd]


def kallisto(inputFastq, outputDir, outputBam, transcriptomeIndex, cpus, inputFastq2=None, size=180, b=200):
    cmd1 = "kallisto quant --bias --pseudobam -b {0} -l {1} -i {2} -o {3} -t {4}".format(b, size, transcriptomeIndex, outputDir, cpus)

    if inputFastq2 is None:
        cmd1 += " --single {0}".format(inputFastq)
    else:
        cmd1 += " {0} {1}".format(inputFastq, inputFastq2)

    cmd1 += " | samtools view -Sb - > {0}".format(outputBam)

    cmd2 = "kallisto h5dump -o {0} {0}/abundance.h5".format(outputDir)

    return [cmd1, cmd2]


def genomeWideCoverage(inputBam, genomeWindows, output):
    cmd = "bedtools coverage -abam -counts -a {0} -b {1} > {2}".format(inputBam, genomeWindows, output)

    return cmd


def calculateFRiP(inputBam, inputBed, output):
    cmd = "cut -f 1,2,3 {0} |".format(inputBed)
    cmd += " bedtools coverage -counts -abam {0} -b - |".format(inputBam)
    cmd += " awk '{{sum+=$4}} END {{print sum}}' > {0}".format(output)

    return cmd


# def macs2CallPeaksATACSeq(treatmentBam, outputDir, sampleName, genome):
#
#     sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}
#
#     cmd = "macs2 callpeak -t {0}".format(treatmentBam)
#     cmd += " --nomodel --extsize 147 -g {0} -n {1} --outdir {2}".format(sizes[genome], sampleName, outputDir)
#
#     return cmd


def bamToBed(inputBam, outputBed):
    cmd = "bedtools bamtobed -i {0} > {1}".format(inputBam, outputBed)

    return cmd


# def centerPeaksOnMotifs(peakFile, genome, windowWidth, motifFile, outputBed):
#
#     cmd = "annotatePeaks.pl {0} {1} -size {2} -center {3} |".format(peakFile, genome, windowWidth, motifFile)
#     cmd += " awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' |"
#     cmd += """ awk -v OFS='\t' -F '\t' '{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }' |"""
#     cmd += " fix_bedfile_genome_boundaries.py {0} | sortBed > {1}".format(genome, outputBed)
#
#     return cmd


def tssAnalysis(inputBam, tssFile, plotsDir, windowWidth, fragmentsize, genome,
                n_clusters, strand_specific, duplicates):
    import os

    cmd = "python {0}/lib/tss_analysis.py {1} {2} {3}".format(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
        inputBam, tssFile, plotsDir
    )
    cmd += " --window-width {0} --fragment-size {1} --genome {2} --n_clusters {3}".format(
        windowWidth, fragmentsize, genome, n_clusters
    )
    if strand_specific:
        cmd += " --strand-specific"
    if duplicates:
        cmd += " --duplicates"

    return cmd


# def plotCorrelations(inputCoverage, plotsDir):
#     import os
#
#     cmd = "python {2}/lib/correlations.py {0} {1}".format(
#         plotsDir,
#         " ".join(["%s"] * len(inputCoverage)) % tuple(inputCoverage),
#         os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
#     )
#
#     return cmd


def getReadType(bamFile, n=10):
    """
    Gets the read type (single, paired) and length of bam file.
    :param bamFile: Bam file to determine read attributes.
    :type bamFile: str
    :param n: Number of lines to read from bam file.
    :type n: int
    :returns: tuple of (readType=string, readLength=int).
    :rtype: tuple
    """
    import subprocess
    from collections import Counter

    try:
        p = subprocess.Popen(['samtools', 'view', bamFile], stdout=subprocess.PIPE)

        # Count paired alignments
        paired = 0
        readLength = Counter()
        while n > 0:
            line = p.stdout.next().split("\t")
            flag = int(line[1])
            readLength[len(line[9])] += 1
            if 1 & flag:  # check decimal flag contains 1 (paired)
                paired += 1
            n -= 1
        p.kill()
    except:
        print("Cannot read provided bam file.")
        raise

    # Get most abundant read readLength
    readLength = sorted(readLength)[-1]

    # If at least half is paired, return True
    if paired > (n / 2):
        return ("PE", readLength)
    else:
        return ("SE", readLength)


def parseBowtieStats(statsFile):
    """
    Parses Bowtie2 stats file, returns series with values.
    :param statsFile: Bowtie2 output file with alignment statistics.
    :type statsFile: str
    """
    import pandas as pd
    import re

    stats = pd.Series(index=["readCount", "unpaired", "unaligned", "unique", "multiple", "alignmentRate"])

    try:
        with open(statsFile) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return stats

    # total reads
    try:
        line = [i for i in range(len(content)) if " reads; of these:" in content[i]][0]
        stats["readCount"] = re.sub("\D.*", "", content[line])
        if 7 > len(content) > 2:
            line = [i for i in range(len(content)) if "were unpaired; of these:" in content[i]][0]
            stats["unpaired"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        else:
            line = [i for i in range(len(content)) if "were paired; of these:" in content[i]][0]
            stats["unpaired"] = stats["readCount"] - int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "aligned 0 times" in content[i]][0]
        stats["unaligned"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if "aligned exactly 1 time" in content[i]][0]
        stats["unique"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if "aligned >1 times" in content[i]][0]
        stats["multiple"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if "overall alignment rate" in content[i]][0]
        stats["alignmentRate"] = re.sub("\%.*", "", content[line]).strip()
    except IndexError:
        pass
    return stats


def parseDuplicateStats(statsFile):
    """
    Parses sambamba markdup output, returns series with values.
    :param statsFile: sambamba output file with duplicate statistics.
    :type statsFile: str
    """
    import pandas as pd
    import re

    series = pd.Series()

    try:
        with open(statsFile) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return series

    try:
        line = [i for i in range(len(content)) if "single ends (among them " in content[i]][0]
        series["single-ends"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
        series["paired-ends"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
        line = [i for i in range(len(content)) if " duplicates, sorting the list...   done in " in content[i]][0]
        series["duplicates"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
    except IndexError:
        pass
    return series


def getPeakNumber(sample):
    """
    Counts number of peaks from a sample's peak file.
    :param sample: A Sample object with the "peaks" attribute.
    :type sampleName: pipelines.Sample
    """
    import subprocess
    import re

    proc = subprocess.Popen(["wc", "-l", sample.peaks], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    sample["peakNumber"] = re.sub("\D.*", "", out)

    return sample
