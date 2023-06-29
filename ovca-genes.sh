#!/bin/env bash
#
#$ -cwd
#
#
#$ -pe smp 8
#$ -l mem_free=12G
#$ -S /bin/bash
#$ -e /scratch/reyesa1/
#$ -o /scratch/reyesa1/

ROOTDIR=`readlink -f .`

bn=`readlink -f . | sed -n 's/.*consensus.*\///p' | sed 's/_con_output//g'`

echo $bn ### save basename (sample and type of consensus)

### prep sites.txt file

cat sites.txt | awk -v OFS="\t" '{print $1, $2-1, $2, $3,$4,$5 }' | bedtools sort > $bn-sites_allcols.bed
cat sites.txt | awk -v OFS="\t" '{print $1, $2-1, $2 }' | bedtools sort > $bn-sites.bed

### do bedtools intersect

bedtools intersect -wb -a $bn-sites_allcols.bed -b ovca_genes_COSMIC.bed > $bn-sites_OvCa_genes.bed

### line count of intermediate/output files

wc -l $bn-sites*
