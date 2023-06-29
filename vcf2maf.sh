#!/bin/env bash
#
#$ -cwd
#
#$ -l mem_free=40G
#$ -pe smp 8
#$ -S /bin/bash
#$ -e /home/yeagern/scratch
#$ -o /home/yeagern/scratch
#export TEMPDIR=/common/yeagern/tmp
#export TMPDIR=/common/yeagern/tmp

while read -r line
do
#line == input file, line by line

perl /common/yeagern/pvacseq_VEP-files/vcf2maf_hg38_ng.pl --input-vcf /common/yeagern/$line".vep.vcf" --output-maf $line".maf" --tumor-id $line --ref-fasta /common/bermanblab/data/public_data/refgenome/HG38_CRCH38/Homo_sapiens_assembly38.fasta --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --inhibit-vep
done < /common/yeagern/PCAWG_vcf2maf_input.txt
