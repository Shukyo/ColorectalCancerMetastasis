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
