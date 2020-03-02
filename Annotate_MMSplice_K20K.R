####################################################################################################################
####################################################################################################################
# Annotate MMSplice outputs.
# Author: Haiying Kong
# Last Modified: 16 November 2019
####################################################################################################################
####################################################################################################################
setwd('/icgc/dkfzlsdf/analysis/B240/kong/Projects/PANCSTRAT')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set paramenters.
raw.dir = 'Lock/WGS/Kipoi/MMSplice/Raw/K20K/'
anno.dir = 'Lock/WGS/Kipoi/MMSplice/Annotated/K20K/'

####################################################################################################################
####################################################################################################################
# Read in GENCODE exon annotation data.
anno = read.table('/icgc/dkfzlsdf/analysis/B240/kong/Reference/GENCODE/Homo_sapiens.GRCh37.75.chr.uniq_exon_nochr_OnlyExon.txt',
                  header=FALSE, sep='\t')
names(anno) = c('Chrom', 'Start', 'End', 'Strand', 'Anno')
anno$Gene_ID = sapply(anno$Anno, function(x)  sub('gene_id ', '', unlist(strsplit(x, '; '))[1])  )
anno$Gene_Symbol = sapply(anno$Anno, function(x)  sub('gene_name ', '', unlist(strsplit(x, '; '))[2])  )
anno $exons = paste0(anno$Chrom, '_', anno$Start, '_', anno$End, ':', anno$Strand)
anno = anno[ ,c('exons', 'Gene_ID', 'Gene_Symbol')]

####################################################################################################################
# Get list of output files from MMSplice.
mmsplice.files = dir(raw.dir, pattern='.txt')

####################################################################################################################
# Clean the result folder.
unlink(anno.dir, recursive=TRUE)
dir.create(anno.dir)

# For each MMSplice output file, annotate and save.
for (mmsplice.file in mmsplice.files)    {
  mmsplice = read.table(paste0(raw.dir, mmsplice.file), header=TRUE, sep='\t')
  mmsplice$ID = sapply(mmsplice$ID, function(x)   paste(unlist(strsplit(x, ':'))[1:2], collapse='_'))

  aster = merge(anno, mmsplice, all.x=FALSE, all.y=TRUE, by='exons')
  names(aster) = gsub('mmsplice_', '', names(aster))

  aster$d_acceptorIntron = aster$alt_acceptorIntron - aster$ref_acceptorIntron
  aster$d_acceptor = aster$alt_acceptor - aster$ref_acceptor
  aster$d_exon = aster$alt_exon - aster$ref_exon
  aster$d_donor = aster$alt_donor - aster$ref_donor
  aster$d_donorIntron = aster$alt_donorIntron - aster$ref_donorIntron

  aster = aster[ ,c('Gene_Symbol', 'Gene_ID', 'dlogitPsi', 'pathogenicity', 'dse',
                    'd_acceptorIntron', 'd_acceptor', 'd_exon', 'd_donor', 'd_donorIntron',
                    'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
		    'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron')]
  aster = aster[order(aster$Gene_Symbol), ]
  new.mmsplice.file = paste0(sub('K20K-', '', unlist(strsplit(mmsplice.file, '_'))[1]), '.txt')
  write.table(aster, paste0(anno.dir, new.mmsplice.file), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }


####################################################################################################################
####################################################################################################################
