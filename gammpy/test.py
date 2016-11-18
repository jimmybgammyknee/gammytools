
while read


def cleanUp(dirtybam, cleanbam):
    """
    Clean up bam file by sorting then removing duplicates
    """

    clean = subprocess.Popen("sambamba sort -p -t {0} {1} | sambamba markdup -p -r -t {0} /dev/stdin {2}".format(threads, dirtybam, cleanbam), shell=TRUE)
    return cleanbam
  
