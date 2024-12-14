library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(tidyverse)
library(Biostrings)
library(purrr)
library(edgeR)
library(ggplot2)
library(WriteXLS)

############## organisation of gene - transcript - exons

#### gene

OrganismDbi::select(x=Mus.musculus,keys='Naca',keytype='SYMBOL',columns='ENTREZID')
OrganismDbi::select(x=Mus.musculus,keys='17938',keytype='ENTREZID',columns='SYMBOL')

genes <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm39.refGene)
genes[genes$gene_id=='17938']


#### exons

exonsBy_gene <- exonsBy(TxDb.Mmusculus.UCSC.mm39.refGene,by = 'gene')
exonsBy_gene$`17938`

#### transcript

txBy_gene <- transcriptsBy(TxDb.Mmusculus.UCSC.mm39.refGene,by='gene')
txBy_gene$`17938`

exonBy_tx <- exonsBy(TxDb.Mmusculus.UCSC.mm39.refGene,by='tx')
exon_tx27443 <- exonBy_tx$`27443`
exon_tx27444 <- exonBy_tx$`27444`
exon_tx27445 <- exonBy_tx$`27445`

exon_tx27443$length <- width(exon_tx27443)
exon_tx27444$length <- width(exon_tx27444)
exon_tx27445$length <- width(exon_tx27445)

############

#transcript <- readDNAStringSet('../0000_database/kallisto_index/Mus_musculus_GRCm39_112/data/Mus_musculus.GRCm39.cdna.all.fa')
transcript <- readDNAStringSet('../00_database/rnaseq/salmon_index/reference_tx/Mus_musculus.GRCm39.cdna.all.fa')

name_tx <- names(transcript)
#map_chr(.x = name_tx,.f = function(x) strsplit(x,split=' c')[[1]][1])

transcript_Naca <- transcript[grepl(pattern = 'Naca ',x=name_tx)]
names(transcript_Naca)
transcript_Naca




################## star count

STARCount_exon <- read.csv('results/STAR_counts_exons.csv')
STARCount_exon %>% 
  filter(GeneID=='17938') 

STARCount_exon_mat <- STARCount_exon %>%
  select(starts_with('SRR'))
STARCount_exon_annotation <- STARCount_exon %>%
  select(GeneID,Length)

dge1 <- DGEList(counts=STARCount_exon_mat,genes =STARCount_exon_annotation )
rownames(dge1$counts) <- paste(dge1$genes$GeneID,dge1$genes$Length,sep='_')
logcpm <- cpm(dge1,prior.count = 1,log = T)

logcpm_selected <- logcpm[substr(rownames(logcpm),1,5)=='17938',]

logcpm_selected_long_exons <- logcpm_selected[4,]
logcpm_selected_long_exons <- data.frame(logcpm_selected_long_exons) %>%
  rownames_to_column(var='Run')

samplelist <- read.csv('data/sample_list.csv')

logcpm_selected_long_exons <- logcpm_selected_long_exons %>%
  left_join(samplelist)

pdf('results/scatterplot_logcpm_long_exon_17938.pdf',width=4,height = 4)
ggplot(data=logcpm_selected_long_exons,mapping = aes(x=Source_name,y=logcpm_selected_long_exons,col=Source_name))+
  geom_jitter(width=0.05,size=0.5)+
  stat_summary(geom="point",fun="mean",size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+
  ylab('Log CPM')+
  xlab('')
dev.off()

WriteXLS(logcpm_selected_long_exons,'results/logcpm_long_exon_17938.xls',row.names = F)





################## salmon count

SalmonCount <- read.csv('results/Salmon_counts_tx.csv')
SalmonCount %>% 
  filter(tx=='ENSMUST00000073868.9'|tx=='ENSMUST00000092048.13') 

SalmonCount <- SalmonCount %>%
  dplyr::select(-X) %>%
  column_to_rownames(var='tx')

dge1 <- DGEList(counts=SalmonCount)
logcpm <- cpm(dge1,prior.count = 1,log = T)

############ ENSMUST00000073868.9

logcpm_selected <- logcpm[rownames(logcpm)=='ENSMUST00000073868.9',]
logcpm_selected <- data.frame(logcpm_selected) %>%
  rownames_to_column(var='Run')
samplelist <- read.csv('data/sample_list.csv')
logcpm_selected <- logcpm_selected %>%
  left_join(samplelist)

pdf('results/scatterplot_logcpm_ENSMUST00000073868.9.pdf',width=4,height = 4)
ggplot(data=logcpm_selected,mapping = aes(x=Source_name,y=logcpm_selected,col=Source_name))+
  geom_jitter(width=0.05,size=0.5)+
  stat_summary(geom="point",fun="mean",size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+
  ylab('Log CPM')+
  xlab('')
dev.off()

WriteXLS(logcpm_selected,'results/logcpm_ENSMUST00000073868.9.xls',row.names = F)


############ ENSMUST00000092048.13

logcpm_selected <- logcpm[rownames(logcpm)=='ENSMUST00000092048.13',]
logcpm_selected <- data.frame(logcpm_selected) %>%
  rownames_to_column(var='Run')
samplelist <- read.csv('data/sample_list.csv')
logcpm_selected <- logcpm_selected %>%
  left_join(samplelist)

pdf('results/scatterplot_logcpm_ENSMUST00000092048.13.pdf',width=4,height = 4)
ggplot(data=logcpm_selected,mapping = aes(x=Source_name,y=logcpm_selected,col=Source_name))+
  geom_jitter(width=0.05,size=0.5)+
  stat_summary(geom="point",fun="mean",size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+
  ylab('Log CPM')+
  xlab('')
dev.off()

WriteXLS(logcpm_selected,'results/logcpm_ENSMUST00000092048.13.xls',row.names = F)







