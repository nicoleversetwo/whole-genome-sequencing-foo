## Dependencies ================================
library(mutSigExtractor)
library(CHORD)
library(BSgenome.Hsapiens.UCSC.hg38)

## Paths ================================

vcf_mutect <- ('~/25236-2315_Mutect2.vcf.gz')
vcf_gridss <- ('~/25236-2315.gripss.somatic.filtered.vcf.gz') ## GRIDSS vcf (SVs), not GRIPSS! (SMNVs)

## Mutation contexts ================================
variants <- list()
contexts <- list()
reference_genome <- BSgenome.Hsapiens.UCSC.hg38
selected_chromosomes <- paste0('chr',c(1:22,'X'))

## SMNVs --------------------------------
## Extract variants
variants$smnv <- mutSigExtractor::variantsFromVcf(
  vcf.file=vcf_mutect, 
  vcf.filter='PASS', 
  keep.chroms=selected_chromosomes, 
  verbose=TRUE
)
#print(variants$smnv)

## Count SNV and indel contexts
contexts$snv <- mutSigExtractor::extractSigsSnv(
  df=variants$smnv, 
  ref.genome=reference_genome, 
  output='contexts'
)

contexts$indel <- mutSigExtractor::extractSigsIndel(
  df=variants$smnv, 
  ref.genome=reference_genome, 
  output='contexts'
)

## SVs --------------------------------
## Extract variants
variants$sv <- mutSigExtractor::variantsFromVcf(
  vcf.file=vcf_gridss, 
  vcf.filter='PASS', 
  keep.chroms=selected_chromosomes,
  vcf.fields=c("CHROM","POS", "REF", "ALT", "FILTER", "ID", "INFO"), 
  verbose=TRUE
)

## Fix chromosome naming from 'chr1' to '1'
variants$sv$chrom <- sub('chr','', variants$sv$chrom) 
variants$sv$alt <- sub('chr','', variants$sv$alt)

##
tmp <- mutSigExtractor::getContextsSv(variants$sv, sv.caller='gridss')
tmp <- mutSigExtractor::extractSigsSv(df=tmp, output='contexts')
contexts$sv <- tmp
rm(tmp)

## Make mutation context matrix for all mutation types --------------------------------
contexts_merged <- do.call(rbind, unname(contexts))
contexts_merged <- t(contexts_merged)

## Predict HRD status ================================
CHORD::chordPredict(contexts_merged)
