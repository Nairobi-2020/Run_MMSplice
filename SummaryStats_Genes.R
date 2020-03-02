####################################################################################################################
####################################################################################################################
# Summary statistics for MMSplice output.
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
summ.dir = 'Lock/WGS/Kipoi/MMSplice/Summary/'

####################################################################################################################
####################################################################################################################
# Clean the result folder.
if (file.exists(summ.dir))  unlink(summ.dir, recursive=TRUE)
dir.create(summ.dir)

####################################################################################################################
# Create large table with counts of SNVs in the genomic regions of our interest.
anno.file = dir(paste0(anno.dir, 'K20K'), pattern='.txt')[1]
aster = read.table(paste0(anno.dir, 'K20K/', anno.file), header=TRUE, sep='\t')[ ,1]
aster = table(aster)
apple = data.frame(Gene_Symbol = names(aster),
                   PID = as.numeric(aster))
names(apple)[2] = sub('.txt', '', anno.file)

for (cohort in c('K20K', 'PCAWG'))   {
  anno.files = dir(paste0(anno.dir, cohort), pattern='.txt')
  if (cohort == 'K20K')  anno.files = anno.files[-1]
  for (anno.file in anno.files)    {
    aster = read.table(paste0(anno.dir, cohort, '/', anno.file), header=TRUE, sep='\t')[ ,1]
    aster = table(aster)
    aster = data.frame(Gene_Symbol = names(aster),
                       PID = as.numeric(aster))
    names(aster)[2] = sub('.txt', '', anno.file)

    apple = merge(apple, aster, by='Gene_Symbol', all.x=TRUE, all.y=TRUE)
    }
  }

for (j in 2:ncol(apple)) apple[is.na(apple[ ,j]), j] = 0

write.table(apple, paste0(summ.dir, 'Counts_SNV_by_GeneTumor.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
# Count SNV in the genomic region of our interest by tumor.
n.snv = apply(apple[ ,-1], 2, sum)
aster = data.frame(PID = names(n.snv),
                   N_SNV = as.numeric(n.snv))

write.table(aster, paste0(summ.dir, 'Counts_SNV_by_Tumor.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
# Count tumors by genes.
n.tumor = apply(apple[ ,-1], 1, function(x) length(which(x!=0)))
aster = data.frame(Gene_Symbol = apple$Gene_Symbol,
                   N_Tumor = n.tumor)
aster = aster[order(-aster$N_Tumor), ]

write.table(aster, paste0(summ.dir, 'Counts_Tumor_by_Gene.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

apple = aster

aster = sapply((aster$N_Tumor)/10, function(x) as.integer(x))
aster = table(aster)
aster = data.frame(Rages = paste0(names(aster), '0s'),
                   N_Genes = as.numeric(aster) )
aster = aster[order(aster$Rages, decreasing=TRUE), ]

write.table(aster, paste0(summ.dir, 'Counts_Genes_by_TumorFreqRange.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

apple = as.data.frame(table(apple$N_Tumor))
names(apple) = c('TumorFreq', 'N_Genes')
apple$TumorFreq = as.numeric(as.character(apple$TumorFreq))
apple = apple[order(-apple$TumorFreq), ]

write.table(apple, paste0(summ.dir, 'Counts_Genes_by_TumorFreq.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')


####################################################################################################################
####################################################################################################################
