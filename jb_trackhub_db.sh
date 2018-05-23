#!/bin/bash -l

bw=$1

echo "track  $(basename $bw _Aligned.sortedByCoord.out.bedgraph.bw)\n
bigDataUrl ../data/${bw}\n
shortLabel $(basename $bw _Aligned.sortedByCoord.out.bedgraph.bw)\n
longLabel $(basename $bw _Aligned.sortedByCoord.out.bedgraph.bw)\n
type bigWig\n
visibility dense\n"
