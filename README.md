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
2. copy __"filter.sh"__to mutect2/filteredVcf  directory

for i in ../rawVcf/*.vcf.gz
do
  name=`basename $i .vcf.gz`
  ./filter.sh $name $i ../../contamination/$name.calculatecontamination.table  ../../ffpe/$name.artifact.pre_adapter_detail_metrics
done


for  i in BXY  DDY  DWY  LDS  WLF  YHD
do
  cd $i
  for j in *.vcf; do convert2annovar.pl --format vcf4old  $j 2>/dev/null |cut -f 1-5  >>$i.raw.bed; done
    sort -V $i.raw.bed|uniq >$i.sorted.uniq.bed
    cd ..
  done

```shell
for k in ../../bamFiles/BXY*.bam
do
    nohup samtools mpileup -f /database/ref/hg19.fa -q 1 -l BXY.sorted.uniq.bed  $i 2>/dev/null >`basename $i .sorted.dedup.rg.recal.bam`.mpileup &
done
```
```shell
for i in */*.mpileup
do
    nohup java -jar /software/VarScan.v2.3.9.jar mpileup2cns $i  --min-coverage  1  --min-reads2  0  --p-value  1  --min-var-freq  0  --strand-filter  0 2>&1 >`dirname $i`/`basename $i .mpileup`.retrieved.vcf &
done
```

./merge_file_from_mpileup2cns.pl *.retrieved.vcf



 for i in */*.somatic_hardfiltered.clean.vcf; do
 nohup gatk HaplotypeCaller -L `dirname $i`/`dirname $i`.sorted.uniq.bed --reference /database/ref/hg19.fa  --input ../bamFiles/`basename $i .somatic_hardfiltered.clean.vcf`.sorted.dedup.rg.recal.bam   --output `dirname $i`/`basename $i .somatic_hardfiltered.clean.vcf`.g.vcf   -ERC GVCF   --output-mode EMIT_ALL_SITES >log.`basename  $i .somatic_hardfiltered.clean.vcf` &
 done

for i in BXY  DDY  DWY  LDS  WLF  YHD
do
   for i in $i/*.g.vcf;do echo " -V "$i;done|xargs  gatk CombineGVCFs -R /database/ref/hg19.fa  -O $i/${i}.g.vcf.gz
done


for i in BXY  DDY  DWY  LDS  WLF  YHD
do
    gatk GenotypeGVCFs -L ${i}/${i}.sorted.uniq.bed  -R /database/ref/hg19.fa -V ${i}/${i}.g.vcf.gz  -O ${i}/${i}.jointlycall.vcf.gz
done


for i in  DDY  DWY  LDS  WLF  YHD
do
  #table_annovar.pl ${i}/${i}.jointlycall.vcf.gz /database/annotation/annovar/humandb/ -buildver hg19 -out ${i}/${i} -remove -protocol refGene,cytoBand -operation g,r -nastring . -vcfinput
  #gatk VariantsToTable -V ${i}/${i}.hg19_multianno.vcf  -F cytoBand -F POS -F REF -F ALT -F Func.refGene -F Gene.refGene -GF AD  -O ${i}/${i}.table
  #awk 'BEGIN{OFS="\t"} {if(NR==1){printf "Chromosome\tPosition\tChange\tGene";for(i=7;i<=NF;i++){split($i,m,".");printf "\t"m[1]}printf "\n"}else{  if ($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/ ) && $NF ~ /,0$/ ){printf "chr"$1"\t"$2"\t"$3">"$4"\t"$6;for(i=7;i<=NF;i++){split($i,m,",");printf "\t"m[2]}printf "\n"}}}' ${i}/${i}.table  >${i}/${i}.mutatedRead.txt
  #awk 'BEGIN{OFS="\t"} {if(NR==1){printf "Chromosome\tPosition\tChange\tGene";for(i=7;i<=NF;i++){split($i,m,".");printf "\t"m[1]}printf "\n"}else{  if ($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/ ) && $NF ~ /,0$/ ){printf "chr"$1"\t"$2"\t"$3">"$4"\t"$6;for(i=7;i<=NF;i++){split($i,m,",");printf "\t"(m[1]+m[2])}printf "\n"}}}' ${i}/${i}.table  >${i}/${i}.coverage.txt
  cp $i/$i.coverage.txt ~/treeomics/src/input/$i/
  cp $i/$i.mutatedRead.txt ~/treeomics/src/input/$i/
done

for i in BXY DDY  DWY  LDS  WLF  YHD
do
awk 'BEGIN{OFS="\t"}{ if(NR==1){for(i=7;i<=NF;i++){split($i,n,".");f[i]=n[1]; printf "" >f[i]".txt"}}  else{ split($1,h,"[pq]"); for(i=7;i<=NF;i++){   split($i,m,",");  if($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/) && length(m)==2  && $NF ~ /,0$/ ){printf "chr"h[1]"\t"($2-1)"\t"$2"\tchr"h[1]":"$2":"$3":"$4"\t"m[1]"\t"m[2]"\t2\t0\t-""\n"  >>f[i]".txt"   }     }  }  }' $i.table
done


for i in *.seg
do
intersectBed -a `basename $i .seg`.txt -b  $i -wb |awk 'BEGIN{OFS="\t"}{if(NR==1){print  "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n"$4,$5,$6,$7,$8,int(2^$14 * 2 +0.5)}else{print $4,$5,$6,$7,$8,int(2^$14 * 2 +0.5)}}' - >`basename $i .seg`.tsv
done

for i in *NOM*.txt
do
awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7,$8,2}' $i |sed "1imutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn">`basename $i .txt`.tsv
done

source activate python27

for i in ../tsv/*.tsv
do
PyClone build_mutations_file --prior total_copy_number --in_file $i --out_file `basename $i .tsv`.yaml
done

for i in bxy ddy dwy lds wlf yhd; do mkdir trace_${i}; done

PyClone run_analysis --config_file config.bxy.yaml --seed 1986




printf "$1:$2:$4:$5\tm[1]\t[m2]"

awk 'BEGIN{OFS="\t"} {if(NR==1){printf "Chromosome\tPosition\tChange\tGene";for(i=7;i<=NF;i++){split($i,m,".");printf "\t"m[1]}printf "\n"}else{  if ($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/ ) && $NF ~ /,0$/ ){printf "chr"$1"\t"$2"\t"$3">"$4"\t"$6;for(i=7;i<=NF;i++){split($i,m,",");printf "\t"m[2]}printf "\n"}}}' BXY.table  >BXY.mutatedRead.txt
awk 'BEGIN{OFS="\t"} {if(NR==1){printf "Chromosome\tPosition\tChange\tGene";for(i=7;i<=NF;i++){split($i,m,".");printf "\t"m[1]}printf "\n"}else{  if ($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/ ) && $NF ~ /,0$/ ){printf "chr"$1"\t"$2"\t"$3">"$4"\t"$6;for(i=7;i<=NF;i++){split($i,m,",");printf "\t"(m[1]+m[2])}printf "\n"}}}' BXY.table  >BXY.coverage.txt
