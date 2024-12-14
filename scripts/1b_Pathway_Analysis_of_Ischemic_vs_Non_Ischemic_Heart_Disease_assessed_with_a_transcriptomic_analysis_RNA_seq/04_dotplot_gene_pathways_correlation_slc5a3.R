library(WebGestaltR)
library(tidyverse)
library(Homo.sapiens)
library(ggplot2)
library(readxl)
library(msigdbr)


#listOrganism()
listGeneSet(organism = 'hsapiens',hostName = "https://www.webgestalt.org/",cache = NULL)

############# KEGG

geneset <- loadGeneSet(organism = 'hsapiens',enrichDatabase='pathway_KEGG')
description <- geneset$geneSetDes
pathway <- geneset$geneSet %>%
  dplyr::select(geneSet,gene) %>%
  left_join(description) %>%
  mutate(geneSet = paste(geneSet,description,sep=' : ')) %>%
  dplyr::select(-description) 

pathway$SYMBOL <- mapIds(x=Homo.sapiens,keys = pathway$gene,keytype = 'ENTREZID',column = 'SYMBOL')

pathway <- pathway %>%
  filter(grepl('hsa04350|hsa04714|hsa04010|hsa04310|04330|04330|04066|04512|04510|04520|04540|04020|04260|04261|04024|04725|04022|04744|04724|04810|05410|04020|04152|05410|05414|04260|04512|00071|00010|00562|04070|00562',geneSet,ignore.case='T'))

table(pathway$geneSet)


################## merge with LOGFC

DGE <- read_excel('../results/Correlation_SLC5A3/DGE.xls')
DGE <- DGE %>%
  dplyr::select(SYMBOL=SYMBOL,logFC,PValue=PValue,FDR=FDR,mlog10PValue)

pathway <- pathway %>%
  left_join(DGE)  

################## merge with GSEA


KEGG_GSEA <- read_excel('../results/Correlation_SLC5A3/KEGG_multigsea/KEGG_GSEA.xls')

KEGG_GSEA <- data.frame(KEGG_GSEA) %>%
  dplyr::select(geneSet=pathway,NES)

pathway <- pathway %>%
  left_join(KEGG_GSEA)  

#########"

pathway <- pathway %>%
  na.omit() %>%
  arrange(logFC) %>%
  mutate(SYMBOL= fct_reorder(SYMBOL,logFC)) %>%
  mutate(geneSet = fct_reorder(geneSet,NES,.desc = T)) %>%
  filter(!is.na(SYMBOL))

pathway <- pathway %>%
  filter(abs(logFC)>0.5) %>%
  filter(FDR<0.01)

pathway$logFC[pathway$logFC>1.5] <- 1.5
pathway$logFC[pathway$logFC<(-1.5)] <- (-1.5)


pdf('../results/Correlation_SLC5A3/KEGG_multigsea/dotplot_gene_pathways_selected.pdf',width=50,height = 7)
 ggplot(data=pathway,mapping = aes(x=SYMBOL,y=geneSet,col=logFC))+
   geom_point(aes(size = mlog10PValue),alpha=0.8)+
   scale_color_gradient2(name=expression(log[2]~(fold-change)),low = "green4", mid = gray(0.7), high = "red3",midpoint=0,limits=c(-1.5,1.5))  +
   theme_bw()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   xlab('')+
   ylab('')+
   labs(size=expression(-log[10]~(P-value)))
dev.off()
 









