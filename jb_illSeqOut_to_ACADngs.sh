#!/bin/bash -l

# Send SequenceOutput directory from demultiplexed Illumina run to ACAD NGS database

# jb_illSeqOut_to_ACADngs.sh [RunName]

name=$1

mkdir -p /data/acad/NGS/"$name"

for file in ./SequenceOutput/*
 do
	fn=$(readlink -f "$file")
	ln -s $fn /data/acad/NGS/"$name"/
done
