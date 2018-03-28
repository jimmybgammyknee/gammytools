#!/bin/bash -l

BAM=$1
THREADS=4
REF=/localscratch/Refs/human/hg19_GRCh37d5/hg19_1000g_hs37d5.fasta
MILLS=/localscratch/Programs/bcbio/genomes/Hsapiens/GRCh37/variation/Mills_and_1000G_gold_standard.indels.vcf.gz
PHASE1=/localscratch/Refs/human/hg19_GRCh37d5/1000G_phase1.indels.hg19.sites.vcf.gz
DBSNP=/localscratch/Programs/bcbio/genomes/Hsapiens/GRCh37/variation/dbsnp-150.vcf.gz
GATK=/apps/software/GATK/3.6-Java-1.8.0_101/GenomeAnalysisTK.jar

picard AddOrReplaceReadGroups I=${BAM}  O=$(basename ${BAM} .bam).rg.bam \
   SO=coordinate RGID=$(basename ${BAM} .bam) RGLB=$(basename ${BAM} .bam) \
   RGPL=Illumina RGPU=NextSeq RGSM=$(basename ${BAM} .bam)

picard MarkDuplicates I=$(basename ${BAM} .bam).rg.bam O=$(basename ${BAM} .bam).markdup.bam \
   M=$(basename ${BAM} .bam).metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

java -jar ${GATK} -T SplitNCigarReads -R ${REF} -I $(basename ${BAM} .bam).markdup.bam \
   -o $(basename ${BAM} .bam).markdup.split.bam -rf ReassignOneMappingQuality \
   -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ${GATK} -T HaplotypeCaller -R ${REF} -I $(basename ${BAM} .bam).markdup.split.bam \
   -o $(basename ${BAM} .bam).markdup.split.vcf.gz --dbsnp ${DBSNP} -dontUseSoftClippedBases \
   -stand_call_conf 30 -stand_emit_conf 30 -nct ${THREADS}

java -jar ${GATK} -T VariantFiltration -R ${REF} -V $(basename ${BAM} .bam).markdup.split.vcf.gz \
   -window 35 -cluster 3 \
   --filterExpression "FS > 30.0 || QD < 2"  -filterName "RNASeqFilters_FS_QD" \
   --filterExpression "QUAL < 3.0 || MQ < 30.0 || DP < 10.0 "  -filterName "LowQualFilter" \
   -o $(basename ${BAM} .bam).markdup.split.filtered.vcf.gz

java -jar ${GATK} -T VariantEval -R ${REF} -nct ${THREADS} \
   -eval $(basename ${BAM} .bam).markdup.split.filtered.vcf.gz -D ${DBSNP} \
   -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants \
   -EV MultiallelicSummary -o $(basename ${BAM} .bam).markdup.split.filtered.vcf.eval.grp

# java -jar ${GATK} -T SelectVariants -R ${REF} -nct ${THREADS} \
#    --variant $(basename ${BAM} .bam).markdup.split.filtered.vcf \
#    -o $(basename ${BAM} .bam).markdup.split.filtered.vcf --selectTypeToExclude INDEL --concordance

java -jar ${GATK} -T ASEReadCounter -R ${REF} -o $(basename ${BAM} .bam).ASE.csv \
   -I $(basename ${BAM} .bam).markdup.bam \
   -sites $(basename ${BAM} .bam).markdup.split.filtered.vcf.gz -U ALLOW_N_CIGAR_READS \
   -minDepth 10 --minMappingQuality 10 --minBaseQuality 2 -drf DuplicateRead
