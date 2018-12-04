#!/bin/bash
for i in ../fastq/*0_NDHE*
do
	name=`expr substr $i 10 8`
	dna_analysis.sh -n  $name -t  16  -s fastq -e reAlign  -1  $i/*1.clean.fq.gz  -2  $i/*2.clean.fq.gz
done
	
