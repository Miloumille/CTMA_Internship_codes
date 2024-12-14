
library(Rsubread)
library(Mus.musculus)
library(tidyverse)
library(openxlsx)
library(edgeR)
library(WriteXLS)


################# 

bam.files <- list.files(path = "data_processed/STAR_mapping/", pattern = ".bam", full.names = TRUE,recursive = T)

#props <- propmapped(files=bam.files)
#props

### export at exon level

fc_exon_level <- featureCounts(bam.files, annot.inbuilt="mm39",isPairedEnd=T,nthreads=8,useMetaFeatures=F)
dim(fc_exon_level$counts)
head(fc_exon_level$counts)
colSums(fc_exon_level$counts)

counts <- fc_exon_level$counts
colnames(counts) <- gsub(colnames(counts) ,pattern='Aligned.sortedByCoord.out.bam',replacement='')
counts_exons <- data.frame(fc_exon_level$annotation,counts)
write.csv(data.frame(counts_exons),'results/STAR_counts_exons.csv',row.names = F)


#################" export at gene level

bam.files <- list.files(path = "data_processed/STAR_mapping/", pattern = ".bam", full.names = TRUE,recursive = T)
fc_gene_level <- featureCounts(bam.files, annot.inbuilt="mm39",isPairedEnd=T,nthreads=8,useMetaFeatures=T)

counts <- fc_gene_level$counts
colnames(counts) <- gsub(colnames(counts) ,pattern='Aligned.sortedByCoord.out.bam',replacement='')
counts_genes <- data.frame(fc_gene_level$annotation %>% dplyr::select(GeneID,Length),counts)
write.csv(data.frame(counts_genes),'results/STAR_counts_genes.csv',row.names = F)




