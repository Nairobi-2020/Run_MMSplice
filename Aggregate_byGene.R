####################################################################################################################
####################################################################################################################
# Aggregate the features of MMSplice output by genes.
# Author: Haiying Kong
# Last Modified: 16 November 2019
####################################################################################################################
####################################################################################################################
setwd('/icgc/dkfzlsdf/analysis/B240/kong/Projects/PANCSTRAT')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set paramenters.
anno.dir = 'Lock/WGS/Kipoi/MMSplice/Annotated/'
aggr.dir = 'Lock/WGS/Kipoi/MMSplice/Aggregated_byGene/Max/'

####################################################################################################################
####################################################################################################################
# Clean the result folder.
if (file.exists(aggr.dir))  unlink(aggr.dir, recursive=TRUE)
dir.create(aggr.dir)
dir.create(paste0(aggr.dir, 'byTumor'))
dir.create(paste0(aggr.dir, 'byFeature'))
dir.create(paste0(aggr.dir, 'byFeature/txt'))
dir.create(paste0(aggr.dir, 'byFeature/forDeep'))

####################################################################################################################
# Read in summary table for tumor frequencies for genes.
summ = read.table('Lock/WGS/Kipoi/MMSplice/Summary/Counts_Tumor_by_Gene.txt', header=TRUE, sep='\t')

####################################################################################################################
####################################################################################################################
# Aggregate features by genes for each tumor.
for (cohort in c('K20K', 'PCAWG'))    {
  anno.files = dir(paste0(anno.dir, cohort), pattern='.txt')

  for (anno.file in anno.files)    {
    aster = read.table(paste0(anno.dir, cohort, '/', anno.file), header=TRUE, sep='\t')[ ,c(1,3:10)]
    apple = data.frame(Gene_Symbol = unique(aster$Gene_Symbol))

    for (j in 2:ncol(aster))    {
      one = aggregate(aster[ ,j], by=list(aster$Gene_Symbol), FUN=function(x) x[(which(abs(x)==max(abs(x))))[1]])
      if (identical(apple$Gene_Symbol, one[ ,1]))
        apple[ ,names(aster)[j]] = one[ ,2]   else
        print('Gene symbols do not match.')
      }

    apple = apple[order(apple$Gene_Symbol), ]
    write.table(apple, paste0(aggr.dir, 'byTumor/', anno.file), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
    }
  }

####################################################################################################################
####################################################################################################################
# Create large table for all tumors for each feature.
tumor.files = dir(paste0(aggr.dir, 'byTumor'), pattern='.txt')
pids = gsub('.txt', '', tumor.files)

features = c('dlogitPsi', 'pathogenicity', 'dse', 'd_acceptorIntron', 'd_acceptor', 'd_exon', 'd_donor', 'd_donorIntron')

for (feature in features)    {
  apple = read.table(paste0(aggr.dir, 'byTumor/', tumor.files[1]), header=TRUE, sep='\t')[ ,c('Gene_Symbol', feature)]
  names(apple)[2] = pids[1]

  for (i in 2:length(tumor.files))    {
    aster = read.table(paste0(aggr.dir, 'byTumor/', tumor.files[i]), header=TRUE, sep='\t')[ ,c('Gene_Symbol', feature)]
    names(aster)[2] = pids[i]
    apple = merge(apple, aster, by='Gene_Symbol', all.x=TRUE, all.y=TRUE)
    }

  for (j in 2:ncol(apple)) apple[is.na(apple[ ,j]), j] = 0

  apple = apple[match(summ$Gene_Symbol, apple$Gene_Symbol), ]
  apple = apple[ ,c('Gene_Symbol', sort(names(apple)[-1]))]
  write.table(apple, paste0(aggr.dir, 'byFeature/txt/', feature, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }


####################################################################################################################
####################################################################################################################
