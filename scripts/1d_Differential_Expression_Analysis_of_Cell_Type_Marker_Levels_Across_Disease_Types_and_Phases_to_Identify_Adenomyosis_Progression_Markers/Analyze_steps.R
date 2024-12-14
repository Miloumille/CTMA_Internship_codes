##Packages
library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)
library(openxlsx)

####################################################################################################

##Raw data
data_multiplex <- read_excel('../data/CD20_CD15_CD117_PanelA_PanelB_JAMBROISE_WACHEUL.xlsx')
data_multiplex$Cells[data_multiplex$Cells=='T cells'] <- 'T_cells'


##Modified data
summary_multiplex_cells <- data_multiplex %>%
  group_by(Image,Type,Cells,Phase,Classification) %>%
  summarise(Result=sum(Result))


##Factor --> ordre données dans graph
summary_multiplex_cells$Type <- factor(summary_multiplex_cells$Type,levels=c('Control','Transition','Adenomyosis'))


##Remove NA si nécessaire
summary_multiplex_cells <- summary_multiplex_cells %>%
  filter(Cells!='NA')



## ?  data_summary  <- data.frame(data_summary,diagnostic=map_chr(.x = data_summary$Classification,.f = function(x) strsplit(x,split='_')[[1]][1]))

####################################################################################################

##GRAPH
pdf('../results/FOLDER_NAME/GRAPH_NAME.pdf', width=6, height=10)
ggplot(summary_multiplex_cells, mapping = aes(
  x = Type,
  y = log2(Result+0.001),
  col = Type)) +
  geom_jitter(width = 0.1, size=0.5, alpha=0.6) +
  stat_summary(geom="point", fun="mean", size=1.5, alpha=0.8) +
  facet_grid(Cells~Phase,scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


####################################################################################################

#### OMIC DATA
summary_wide <- summary_multiplex_cells %>%
  pivot_wider(names_from = Cells,
              values_from = Result)

##Ajout d'un "ID" à chaque data
summary_wide <- data.frame(ID=paste0('S',seq(1,dim(summary_wide)[1])),summary_wide)

## Retrait des "coldata" ((-Image,-Type,-Phase,-Classification)) pour ne garder que les data omic (=Cells + Result)
omic_data <- summary_wide %>%
  ungroup() %>%
  dplyr::select(-Image,-Type,-Phase,-Classification) %>%
  column_to_rownames(var='ID')

omic_data <- t(omic_data)

##Transformation en log2
omic_data <- log2(omic_data+0.001)

## coldata
coldata <- summary_wide %>%
  ungroup() %>%
  dplyr::select(ID,Image,Type,Phase) 



####################################################################################################

##Design 
######## unique sample per patients [Ici : Type]
mydesign <- model.matrix(~Type, data = coldata)

fit <- lmFit(omic_data, design = mydesign)

######## multiple sample per patients [Ici : Classification]
mydesign <- model.matrix(~Classification, data = coldata)

corfit <- duplicateCorrelation(omic_data, design = mydesign, block = coldata$Image)  ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS
fit <- lmFit(omic_data, design = mydesign, block = coldata$Image, correlation = corfit$consensus)   ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS


####################################################################################################

##Comparaison Adenomyosis_vs_control
cm <- makeContrasts(Adenomyosis_vs_control = TypeAdenomyosis, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


###Comparaison Transition_vs_control
cm <- makeContrasts(Transition_vs_control = TypeTransition, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

#### si on veut faire la différence entre Adenomyosis et Transition 
#### on remplace "TypeAdenomyosis" par "TypeAdenomyosis-TypeTransition"


####################################################################################################


####Resultats avec pvalue
result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

####################################################################################################

###Export results
write.xlsx(result,'../results/FOLDER_NAME/FILE_NAME.xlsx')


