library(tximport)
library(Biostrings)
library(edgeR)
library(ggplot2)
library(WriteXLS)
library(purrr)

transcriptome <- readDNAStringSet('../00_database/rnaseq/salmon_index/reference_tx/Mus_musculus.GRCm39.cdna.all.fa')
transcriptnames <- names(transcriptome)
transcriptnames <- unlist(lapply(strsplit(transcriptnames,split = ' cdna'), function(l) l[[1]]))
genenames <- names(transcriptome) 
genenames <- unlist(lapply(strsplit(genenames,split = 'gene:'), function(l) l[[2]]))
genenames <- unlist(lapply(strsplit(genenames,split = ' gene_'), function(l) l[[1]]))
tx2gene <- data.frame(transcriptnames,genenames)

files <- list.files('data_processed/salmon_quant/',full.names = T,pattern = 'quant.sf',recursive = T)



############### estimate counts at transcript levels 

txi <- tximport(files, type = "salmon", tx2gene = tx2gene,txOut=T)
names(txi)
counts <- txi$counts
colnames(counts) <- basename(list.dirs('data_processed/salmon_quant/',recursive = F)) 
counts <- data.frame(counts) %>% rownames_to_column(var='tx')
write.csv(counts,'results/Salmon_counts_tx.csv')



txi <- tximport(files, type = "salmon", tx2gene = tx2gene,txOut=F)
names(txi)
counts <- txi$counts
colnames(counts) <- basename(list.dirs('data_processed/salmon_quant/',recursive = F)) |> substr(1,10)
counts <- data.frame(counts) %>% rownames_to_column(var='ENSEMBL')

write.csv(counts,'results/Salmon_counts_gene.csv')


