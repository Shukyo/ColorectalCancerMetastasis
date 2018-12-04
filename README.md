# ColorectalCancerMetastasis
---
Scripts for processing WES data from Colorectal Cancer Metastasis

## Prerequisite

## Run the scripts
### Alignment and Bam processing
1. copy __"align.sh"__ to working directory
```shell
for name in `ls -d *[0]`
do
  #dna_analysis.sh -n  $name  -t  16  -s fastq -e reCal  -1  $fq1  -2  $fq2
  dna_analysis.sh -n  $name  -t  10  -s reCal -e varCall  -i ${name}/${name}.sorted.dedup.rg.bam
done
```
