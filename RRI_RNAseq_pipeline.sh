#!/bin/bash -l

## Robinson Research Institute RNAseq pipeline

## Jimmy Breen (jimmymbreen@gmail.com)

## Requires:
##	HISAT2
##	FastQC
##	AdapterRemoval
##	featureCounts
##	StringTie
##	Cufflinks (for downstream assembly)

### DEVELOPMENT 
## -Convert to ruffus python pipeline
## -Config file in YAML format
## -Phoenix usage
## -Stable log files

# Module loading
. /opt/shared/Modules/3.2.7/init/bash
module load cufflinks/2.1.1
module load gnu/4.9.2
module load java/java-1.7.09
module load fastQC/0.11.2
#module load htseq/0.6.1p1 ## DEPRECIATED
module load parallel

## Directory setup
base=$(pwd)
DATA=$base/Data
ASS=$base/4_stringtie
QUANT1=$base/3_featureCounts
QUANT2=$base/3_salmon
BAMS=$base/2_Hisat2_merged
TRIM=$base/1_AdapterRemoval
QC=$base/0_FastQC

## Custom parameters and programs
## - Turn this into YAML config file for python pipeline
ar=/localscratch/Programs/adapterremoval2-20160217-git~cf3b668/build/AdapterRemoval
bamba=$HOME/bin/sambamba
featcounts=/localscratch/Programs/subread-1.5.1-Linux-x86_64/bin/featureCounts
stringtie=/localscratch/Programs/stringtie-1.3.0.Linux_x86_64/stringtie
salmon=/localscratch/Programs/Salmon-0.7.2_linux_x86_64/bin/salmon
threads=32

# HISAT2 + transcript References
hg38=/localscratch/Refs/human/hg38_hg20_GRCh38p5/hisat/grch38_tran/genome_tran
#gff=/localscratch/Refs/human/hg38_hg20_GRCh38p5/gencode.v24.chr_patch_hapl_scaff.annotation.gff3
mm10=/localscratch/Refs/mouse/grcm38_snp_tran/genome_snp_tran
gff=/localscratch/Refs/mouse/gencode.vM11.chr_patch_hapl_scaff.annotation.gtf
hg38_trans=/localscratch/Refs/human/hg38_hg20_GRCh38p5/gencode.v24.lncRNA_transcripts_plus_ERCC.fa.gz
hg38_map=
#mm10_map=/localscratch/Refs/mouse/mm38.genenameMap.map
mm10_trans=/localscratch/Refs/mouse/gencode.vM11.transcripts.fa.gz
z10=/localscratch/Refs/Danio_rerio/GRCz10/Danio_rerio.GRCz10.dna_sm.toplevel
z10_gff=/localscratch/Refs/Danio_rerio/GRCz10/Danio_rerio.GRCz10.cdna.all.fa.gz

# QC using FastQC
mkdir -p $QC
for i in $DATA/*.fastq.gz
 do
	fastqc -t $threads -k 9 $i -o $QC/
done

# AdapterTrimming with AdapterRemoval
mkdir -p $TRIM
for FQGZ in $DATA/*R1_001*.fastq.gz
 do
	$ar --file1 $FQGZ --file2 ${FQGZ/R1_001/R2_001} \
		--output1 $TRIM/$(basename $FQGZ _R1_001.fastq.gz)_trim1.fastq.gz \
		--output2 $TRIM/$(basename $FQGZ _R1_001.fastq.gz)_trim2.fastq.gz \
		--threads $threads --gzip \
		--trimqualities --trimns --minquality 10 --minlength 25
	parallel -j"$threads" "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: *_trim*.fastq.gz \
		> $TRIM/project_AdapterRemoval.txt
done

# Alignment with HISAT2
mkdir -p $BAMS
for FQGZ in $TRIM/*_trim1*.fastq.gz
 do 
	hisat2 -p $threads -x $mm10 \
		-1 ${FQGZ} \
		-2 ${FQGZ/trim1/trim2} | \
		samtools view -bhS -q30 - 1> $BAMS/$(basename $FQGZ _trim1.fastq.gz).bam \
			2>> $BAMS/project_alignment.stats
	$bamba sort -p -t $threads -n \
		-o $BAMS/$(basename $FQGZ _trim1.fastq.gz)_hisat2_mm10.sorted.bam \
		$BAMS/$(basename $FQGZ _trim1.fastq.gz)_hisat2_mm10.bam  
done

#Quantification with featureCounts
mkdir -p $QUANT1

samp_list=`find $BAMS -name "*.sorted.bam" | tr '\n' ' '`

featureCounts -Q 10 -s 1 -T $threads -a $gff -o $QUANT1/project_genes.out $samp_list
cut -f1,7- $QUANT1/project_genes.out | sed 1d > $QUANT1/project_genes.txt

#Quantification with salmon
mkdir -p $QUANT2

$salmon index -p $threads -i $QUANT2/transcripts_salmon_idx -t $mm10_trans --gencode --perfectHash

for FQGZ in $TRIM/*_trim1*.fastq.gz
 do
	$salmon quant -l A -i $QUANT2/transcripts_salmon_idx -p $threads -1 $FQGZ -2 ${FQGZ/trim1/trim2} --numBootstraps 100 -o $QUANT2/$(basename $FQGZ _trim1.fastq.gz).salmonOut
done

# Assembly
mkdir -p $ASS
for bam in $BAMS/*.sorted.bam 
 do
	$stringtie -p $threads -m 30 -G $gff \
          -o $ASS/$(basename $bam _hisat2_mm10.sorted.bam)_assembly.gtf \
          -A $ASS/$(basename $bam _hisat2_mm10.sorted.bam)_gene_abund.tab \
          -C $ASS/$(basename $bam _hisat2_mm10.sorted.bam)_cov_refs.gtf $bam
done

$stringtie --merge -G $gff -o $ASS/merged.gtf `find $ASS -name "*assembly.gtf"`
awk '{if($3=="transcript")print}' $ASS/merged.gtf > $ASS/merged.transcript.gtf

# from here we cuffcompare with other sets  

