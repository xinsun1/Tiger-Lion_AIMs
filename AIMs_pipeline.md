# Pipeline for AIMs identification

### 1. Genotype calling
#### 1.1 samtools 



### 1. Genotype calling
#### 1.1 samtools 


### 2. Create masks
#### 2.1 repeat elements 



#### 2.2 mappability mask 

``` bash
#!/bin/bash

WDIR=/projects/mjolnir1/people/gnr216/3-panthera_aims/0-ref

cd $WDIR

# careful about memory consumption
# /projects/mjolnir1/people/gnr216/a-software/genmap/genmap-build/bin/genmap index -F ./felcat9.fasta -I ./genmap_index


/projects/mjolnir1/people/gnr216/a-software/genmap/genmap-build/bin/genmap map -K 100 -E 2 -I ./genmap_index -O ./genmap_out -t -w -bg -T 20
```
filter by threshold < 1



#### 2.3 numt mask, (blastn)
run blastn against tiger, lion and cave lion mitochondrial sequences
merge the result


|Species|MT used for blastb|
|---|---|
|Cave lion|KX258452.1|
|Tiger|KP202268.1|
|Lion|KP202262.1|



``` bash
# 1. create blastdb
makeblastdb -in ../0-ref/felcat9.fasta -dbtype nucl -out blastdb_felcat9

# 2. run blastn
blastn -db ./blastdb_felcat9/blastdb_felcat9 -query mt.all.fa -out out.blastn.W16E1e-4.mt_all -word_size 16 -evalue 0.0001 -outfmt 7

# 3. clean result
grep -v -e "chrMT" -e "^#" out.blastn.W16E1e-4.mt_all | awk '{print $2,$9,$10}' | awk '$2>$3 {print $1,$3,$2;next}{print $0}' OFS='\t'
 > out.mt_all.hits.bed

sort -k1,1V -k2,2g out.mt_all.hits.bed | awk '{print $1,$2-200,$3+200}' OFS='\t' | bedtools merge > out.mt_all.hits.f100_merge.bed

```


#### 2.4 ancient DNA mask
maybe try to remove transitions


#### from Jazmin
These are the parameters we used for SNP calling in ANGSD

angsd -bam LINE.bamlist -minq 30 -minmapq 30 -GL 2 -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doglf 2 -minind 5 -minIndDepth 3 -out LINE.allSamples.30pct

After this, our only filter was that for SNPs that were too close (<1000bp) from each other we randomly selected one of them, because we expected a lot of SNPs to be very close in the data we had. 

We decided not to apply a MAF filter, because we had less than a thousand SNPs with MAF below 0.01
if you are using ANGSD for SNP calling, the advantage is that it normally writes a table  with info on the missingness and maf of each SNP, so you can use it to decide the filters later if you end up with too many SNPs



#### 2.5 INDEL mask
avoid regions close to an INDEL



#### 2.6 mask all

``` bash
cat ../0-ref/felcat.rep.trf.bedpos.bed ../5-numt/out.mt_all.hits.f100_merge.bed | sort -k1,1V -k2,2g | bedtools merge > felcat.rep.trf.numt_f100.bedpos.bed

bedtools intersect -a gl_tv_maf05_mis50.bed -b felcat.rep.trf.numt_f100.bedpos.bed -v > gl_tv_maf05_mis50.rep.trf.numt_f100.bed
```

### 3. pop gen
Idea: only keep SNPs that retains the tiger/lion/cave lion structure.
PCA, Fst, AIMs Info (ref: Informativeness of Genetic Markers for Inference of Ancestry)

#### 3.1 pca


#### 3.2 

### 4. check performance

#### 1.1 samtools 


