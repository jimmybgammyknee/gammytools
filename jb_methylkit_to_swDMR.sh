#!/bin/bash -l

CG=$1
CHG=$2
CHH=$3

grep -v "^chrBase" $CG | awk '{cs=($5 * $6); \
	ts=($5 * $7); \
	if ($4 == "F") str= "+"; \
	else if ($4 == "R") str= "-"; \
	print $2"\t"$3"\t"$3+1"\t""CG""\t"$5"\t"int(cs)"\t"int(ts)"\t"str;
	}' > ${CG}.cg.swdmr.tsv

grep -v "^chrBase" $CHG | awk '{cs=($5 * $6); \
        ts=($5 * $7); \
        if ($4 == "F") str= "+"; \
        else if ($4 == "R") str= "-"; \
        print $2"\t"$3"\t"$3+2"\t""CHG""\t"$5"\t"int(cs)"\t"int(ts)"\t"str;
        }' > ${CG}.chg.swdmr.tsv 

grep -v "^chrBase" $CHH | awk '{cs=($5 * $6); \
        ts=($5 * $7); \
        if ($4 == "F") str= "+"; \
        else if ($4 == "R") str= "-"; \
        print $2"\t"$3"\t"$3+2"\t""CHH""\t"$5"\t"int(cs)"\t"int(ts)"\t"str;
        }' > ${CG}.chh.swdmr.tsv
