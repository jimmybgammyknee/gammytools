#!/bin/bash -l

. /opt/shared/Modules/3.2.7/init/bash
module load java/java-jdk-1.7.051
module load gatk/3.1-1

bam=$1
gatk=/opt/shared/gatk/3.1-1-g07a4bf8/GenomeAnalysisTK.jar
ref=/localscratch/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta

java -jar $gatk -R $ref -T UnifiedGenotyper -I $bam -o ${bam}.gatk_GT.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 -rf BadCigar
