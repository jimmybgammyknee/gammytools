#!/bin/bash -l

miseq=$1

FCID=$(sed -n 's:.*<Flowcell>\(.*\)</Flowcell>.*:\1:p' RunInfo.xml)

echo "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject" > SampleSheet_JB.csv

sed -n '/^Sample_ID/,$p' $miseq | \
	awk -v var="$FCID" '{print var",""1"","$1","$2","$6","$5",""N"","",""JB"","$7}' | \
	sed '/Sample_Plate,Sample_Well,I7_Index_ID/d' >> SampleSheet_JB.csv  

