#!/bin/bash

if [ ! $1 ]||[ ! $2 ]||[ ! $3 ]||[ ! $4 ];then
   echo "Usage:"
   echo "     "$0 "prefix inputVcf contamination/calculatecontamination.table ffpe/artifact.pre_adapter_detail_metrics "
   echo
   exit 0
fi

gatk FilterMutectCalls \
        -V $2 \
        --contamination-table   $3  \
        -O $1.somatic_oncefiltered.vcf.gz

gatk FilterByOrientationBias \
        -V $1.somatic_oncefiltered.vcf.gz \
        -P $4  \
        -O $1.somatic_twicefiltered.vcf.gz


zcat $1.somatic_twicefiltered.vcf.gz|awk '/^#/ || $7~/PASS/' - >$1.somatic_twicefiltered.clean.vcf




## Depth >=3 ; Allele Frquency >= 10% ; Variants on both strand
awk -F "\t" '{split($10,f,":");split(f[5],m,",");split(f[6],n,",");if(/^#/ || (f[4]>=3 && f[3] >=0.1 && m[2] > 0 && n[2] > 0 )){print }}'  $1.somatic_twicefiltered.clean.vcf >$1.somatic_hardfiltered.clean.vcf
