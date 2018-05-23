#!/bin/bash -l

# Adjusts BAM read-ends by Tn5 offsets (+4 for forward, and -5 for reverse)
tn5_adjust='        # tn5_adjust is an awk function
BEGIN {OFS = FS} {  # use the set the original files field separator
  if ($2 == 0) {    # if the read is on the positive strand
    $4 = $4 + 4     # shift the alignment upstream 4bp
  } else {          
    $4 = $4 - 5     # else shift it downstream 5bp
  }
  if ($4 < 1) {     # if we hit the beginning of the chromosome do nothing
      $4 = 1
  }
  print $0          # return all rows
}'

SAMPLE=$1                                    # The input BAM   
THREADS=$2                                   # Number of threads for samtools

samtools view -H $1 > ${SAMPLE}.tmp          # Store the original header

samtools view $1 | awk -F $'\t' "$tn5_adjust" >> ${SAMPLE}.tmp  # Run the awk function appending to the header file

samtools view -hb -@ $THREADS \              # Revert SAM back to BAM
	-o ${SAMPLE%.bam}_tn5Adj.bam ${SAMPLE}.tmp   # The output file

rm ${SAMPLE}.tmp                             # remove the tmp file
exit 0
