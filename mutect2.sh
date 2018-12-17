#!/bin/bash

gatk Mutect2  \
	-R /database/ref/hg19.fa \
	-I $1 \
	-I $2 \
	-tumor `basename $1 .sorted.dedup.rg.recal.bam` \
	-normal `basename $2 .sorted.dedup.rg.recal.bam` \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	-O `basename $1 .sorted.dedup.rg.recal.bam`.mutect2.vcf.gz



gatk CreateSomaticPanelOfNormals \
  --min-sample-count 1 \
  -vcfs BXYNOM58.normal.vcf.gz  \
  -vcfs  DDYNOM01.normal.vcf.gz  \
  -vcfs  DWYNOM69.normal.vcf.gz  \
  -vcfs  LDSNOM17.normal.vcf.gz  \
  -vcfs  WLFNOM44.normal.vcf.gz \
  -vcfs   YHDNOM36.normal.vcf.gz \
  -O 6samples_pon.vcf.gz


gatk -java-options "-Xmx2g" Mutect2 \
  -R /database/ref/hg19.fa \
  -I tumor.bam \
  -I normal.bam \
  -tumor HCC1143_tumor \
  -normal HCC1143_normal \
  -pon resources/chr17_pon.vcf.gz \
  --germline-resource resources/chr17_af-only-gnomad_grch38.vcf.gz \
  --af-of-alleles-not-in-resource 0.0000025 \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  -L chr17plus.interval_list \
  -O 1_somatic_m2.vcf.gz \
  -bamout 2_tumor_normal_m2.bam
