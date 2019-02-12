# ColorectalCancerMetastasis
---
Scripts for processing WES data from Colorectal Cancer Metastasis

## Prerequisite


## Alignment and Bam processing
### Alignment
copy __"align.sh"__ to working directory
```shell
for name in `ls -d *[0]`
do
  dna_analysis.sh -n  $name  -t  16  -s fastq -e reCal  -1  $fq1  -2  $fq2
done
```
### Filters
copy __"filter.sh"__ to mutect2/filteredVcf  directory

```shell
for i in ../rawVcf/*.vcf.gz
do
  name=`basename $i .vcf.gz`
  ./filter.sh $name $i ../../contamination/$name.calculatecontamination.table  ../../ffpe/$name.artifact.pre_adapter_detail_metrics
done
```

3. Combine variants for outputing bed file

## Retrieve the selected variants
### Get the bed file for each patient
```shell
for  i in BXY  DDY  DWY  LDS  WLF  YHD
do
  cd $i
  for j in *.vcf
  do
       convert2annovar.pl --format vcf4old  $j 2>/dev/null |cut -f 1-5  >>$i.raw.bed
  done
  sort -V $i.raw.bed|uniq >$i.sorted.uniq.bed
  cd ..
  done
```

### Call gVcfs
```shell
 for i in */*.somatic_hardfiltered.clean.vcf
 do
 nohup gatk HaplotypeCaller -L `dirname $i`/`dirname $i`.sorted.uniq.bed --reference /database/ref/hg19.fa  --input ../bamFiles/`basename $i .somatic_hardfiltered.clean.vcf`.sorted.dedup.rg.recal.bam   --output `dirname $i`/`basename $i .somatic_hardfiltered.clean.vcf`.g.vcf   -ERC GVCF   --output-mode EMIT_ALL_SITES >log.`basename  $i .somatic_hardfiltered.clean.vcf` &
 done
```
### Process the output gVcf
```shell
for i in BXY  DDY  DWY  LDS  WLF  YHD
do
   for i in $i/*.g.vcf;do echo " -V "$i;done|xargs  gatk CombineGVCFs -R /database/ref/hg19.fa  -O $i/${i}.g.vcf.gz
   gatk GenotypeGVCFs -L ${i}/${i}.sorted.uniq.bed  -R /database/ref/hg19.fa -V ${i}/${i}.g.vcf.gz  -O ${i}/${i}.jointlycall.vcf.gz
done
```

### annotation + treeomics Input file prepare

```shell
for i in  BXY DDY  DWY  LDS  WLF  YHD
do
  table_annovar.pl ${i}/${i}.jointlycall.vcf.gz /database/annotation/annovar/humandb/ -buildver hg19 -out ${i}/${i} -remove -protocol refGene,cytoBand -operation g,r -nastring . -vcfinput
  gatk VariantsToTable -V ${i}/${i}.hg19_multianno.vcf  -F cytoBand -F POS -F REF -F ALT -F Func.refGene -F Gene.refGene -GF AD  -O ${i}/${i}.table
  awk 'BEGIN{OFS="\t"} {if(NR==1){printf "Chromosome\tPosition\tChange\tGene";for(i=7;i<=NF;i++){split($i,m,".");printf "\t"m[1]}printf "\n"}else{  if ($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/ ) && $NF ~ /,0$/ ){printf "chr"$1"\t"$2"\t"$3">"$4"\t"$6;for(i=7;i<=NF;i++){split($i,m,",");printf "\t"m[2]}printf "\n"}}}' ${i}/${i}.table  >${i}/${i}.mutatedRead.txt
  awk 'BEGIN{OFS="\t"} {if(NR==1){printf "Chromosome\tPosition\tChange\tGene";for(i=7;i<=NF;i++){split($i,m,".");printf "\t"m[1]}printf "\n"}else{  if ($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/ ) && $NF ~ /,0$/ ){printf "chr"$1"\t"$2"\t"$3">"$4"\t"$6;for(i=7;i<=NF;i++){split($i,m,",");printf "\t"(m[1]+m[2])}printf "\n"}}}' ${i}/${i}.table  >${i}/${i}.coverage.txt
  cp $i/$i.coverage.txt ~/treeomics/src/input/$i/
  cp $i/$i.mutatedRead.txt ~/treeomics/src/input/$i/
done
```
