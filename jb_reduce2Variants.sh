#!/bin/bash -l

# Take "per-bp" gz compressed vcfs (in one directory) that has been split by chromosome 
# and then reduce into one file that only contains variants 

module load vcftools

dir=$1

for i in $dir/*.vcf.gz
 do
	samplename=$(basename $i .vcf.gz)
	tabix -p vcf $i
done 

vcf-concat $dir/*.vcf.gz | awk '$5!="." {print $0}' | bgzip -c > $samplename_all.vcf.gz