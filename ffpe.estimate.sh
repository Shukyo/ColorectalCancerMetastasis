#!/bin/bash
for i in *1.sorted.dedup.realign.recal.bam
do


        # collect metrics on sequence context artifacts
        gatk CollectSequencingArtifactMetrics \
                -I $i \
                -O `expr substr $i 1 8`.artifact \
                -R /database/ref/hg19.fa


        #perform orientation bias filtering
        gatk FilterByOrientationBias \
                -A G/T \
                -A C/T \
                -V 9_somatic_oncefiltered.vcf.gz \
                -P $name.artifact.pre_adapter_detail_metrics \\
                -O $name.somatic_twicefiltered.vcf.gz


done
