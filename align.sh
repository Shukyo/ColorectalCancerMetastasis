#!/bin/bash

name="$1"
fq1="$2"
fq2="$3"

dna_analysis.sh -n  $name  -t  16  -s fastq -e reCal  -1  $fq1  -2  $fq2
