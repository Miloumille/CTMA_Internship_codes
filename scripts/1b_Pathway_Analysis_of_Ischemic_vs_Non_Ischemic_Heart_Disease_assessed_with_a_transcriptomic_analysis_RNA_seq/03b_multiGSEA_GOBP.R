
library(multiGSEA)
library(magrittr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(WriteXLS)
library(Homo.sapiens)


############## put old kegg results in old directory

############### 


## prepare pathway

library(WebGestaltR)
library(tidyverse)
listOrganism()
listGeneSet(organism = 'hsapiens',hostName = "https://www.webgestalt.org/",cache = NULL)
geneset <- loadGeneSet(organism = 'hsapiens',enrichDatabase='geneontology_Biological_Process_noRedundant')

description <- geneset$geneSetDes

pathway <- geneset$geneSet %>%
  dplyr::select(geneSet,gene) %>%
  left_join(description) %>%
  mutate(geneSet = paste(geneSet,description,sep=' : ')) %>%
  dplyr::select(-description) %>%
  unstack(gene ~ geneSet)

pathways <- list("transcriptome" = pathway)

list_results <- list.files('../results',pattern = 'DGE',recursive = T,full.names = T)

for(i in 1:length(list_results))
{
  dir.create(gsub(list_results[i],pattern='DGE.xlsx',replacement='GOBP_multigsea'))
  dge_result <- read_excel(list_results[i])
  dge_result <- data.frame(ENTREZID=mapIds(x=Homo.sapiens,keys = dge_result$SYMBOL,keytype = 'SYMBOL',column = 'ENTREZID'),dge_result)
  
  dge_result <- dge_result %>% filter(!is.na(logFC)) %>%
    filter(!is.na(PValue))
  
  omics_data <- initOmicsDataStructure(layer = c("transcriptome"))
  omics_data$transcriptome <- rankFeatures(logFC=dge_result$logFC,pvalues=dge_result$PValue)
  
  names(omics_data$transcriptome) <- dge_result$ENTREZID
  
  enrichment_scores <- multiGSEA(pathways, omics_data)
  enrichment_scores <- data.frame(enrichment_scores$transcriptome)
  enrichment_scores <- enrichment_scores %>% arrange(pval) %>%
    dplyr::select(-log2err,-leadingEdge)
  
  WriteXLS(enrichment_scores,gsub(list_results[i],pattern='DGE',replacement='GOBP_multigsea/GOBP_GSEA'))
  
}






