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

### PyClone input preparation
> Should run the initial two steps of treeomics input preparation to get table files

#### Split the merged table file into individual files
> ONLY exonic and splicing variants reatined, sites on chrM or with alterantive reads in NORMAL sample were removed

```shell
for i in BXY DDY  DWY  LDS  WLF  YHD
do
awk 'BEGIN{OFS="\t"}{ if(NR==1){for(i=7;i<=NF;i++){split($i,n,".");f[i]=n[1]; printf "" >f[i]".txt"}}  else{ split($1,h,"[pq]"); for(i=7;i<=NF;i++){   split($i,m,",");  if($1 != "." && ($5 ~/^exonic$/ || $5 ~/^splicing$/) && length(m)==2  && $NF ~ /,0$/ ){printf "chr"h[1]"\t"($2-1)"\t"$2"\tchr"h[1]":"$2":"$3":"$4"\t"m[1]"\t"m[2]"\t2\t0\t-""\n"  >>f[i]".txt"   }     }  }  }' $i.table
done
```

#### Get the copy number info
```shell
for i in *.seg
do
  # the copy number MUST be integrate
  intersectBed -a `basename $i .seg`.txt -b  $i -wb |awk 'BEGIN{OFS="\t"}{if(NR==1){print  "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n"$4,$5,$6,$7,$8,int(2^$14 * 2 +0.5)}else{print $4,$5,$6,$7,$8,int(2^$14 * 2 +0.5)}}' - >`basename $i .seg`.tsv
done
  # Missing normal files in the above calculation, create tsv files for normal samples
for i in *NOM*.txt
do
    awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7,$8,2}' $i |sed "1imutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn">`basename $i .txt`.tsv
done
```

#### Run the Pyclone
>Create the config file maunally according the demo

```shell
  # activate the vitrual python enviorment
  source activate python27
  # build the yaml files
  for i in ../tsv/*.tsv
  do
      PyClone build_mutations_file --prior total_copy_number --in_file $i --out_file `basename $i .tsv`.yaml
  done
  # run the main program of PyClone
  for i in bxy ddy dwy lds wlf yhd
  do
      mkdir trace_${i}
      PyClone run_analysis --config_file config.${i}.yaml --outfile trace_${i}   --seed 1986
      PyClone build_table  --config_file config.${i}.yaml --outfile ${i}.table --table_type loci
  done
```

#### Create gene-based tsv files
> Use the tsv files created above;
Use homemake script intersect

```shell
  for i in *.tsv
  do
    name=`basename $i .tsv`
    # Prepare the annotation input files
    tail -n +2 $i |awk -F "[:\t]" 'BEGIN{OFS="\t"}{print $1,$2-length($3)+1,$2,$3,$4,$0}' - >$name.avinput
    # Annotation
    annotate_variation.pl -buildver hg19 -outfile $name $name.avinput /database/annotation/annovar/humandb
    # Remove synonymous variants
    awk '$2~/^synonymous/' $name.exonic_variant_function >$name.tmp
    intersect  -c 8 9  -r  $name.variant_function -f $name.variant_function $name.tmp -o 1 >$name.txt
    rm $name.tmp
    # Create tsv files (should create tsv folder in case overlap original tsv)
    mkdir tsv
    sort -k2  $name.txt |awk 'BEGIN{OFS="\t"}{if($2~/\(/){sub(/\(\S+\)/,"",$2)};m[$2]=m[$2]+$9;n[$2]=n[$2]+$10;h[$2]=h[$2]+$13;a[$2]=a[$2]+1}END{for(k in m){print  k,m[k],n[k],2,0,int(h[k]/a[k])}}' -  >tsv/$name.tsv
done
```
