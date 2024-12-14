##Packages
library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)
library(openxlsx)

####################################################################################################
##Comparaison lesions_AD_versus_transition

##Raw data
cells_multiplex <- read_excel('../data/cell_names_2.xlsx')

##Modified data
AD <- cells_multiplex %>%
  group_by(Phase,Image,Classification,Cells_group) %>%
  summarise(Result=sum(Result)) %>% 
  filter(Classification %in% c("AD_ectopic_adenomyosis", "Transition_ectopic")) %>%
  filter(Cells_group!='NA', Cells_group!='Neutrophils')


##Factor --> ordre données dans graph
AD$Classification <- factor(AD$Classification, levels = c('Transition_ectopic', 'AD_ectopic_adenomyosis'))

####################################################################################################

##GRAPH
#PC PORTABLE

AD$Cells_subgroup <- factor(AD$Cells_subgroup,
                            levels=c('T cells',
                                     'Cytotoxic T cells','Helper 1 T cells',
                                     'Helper 2 T cells','NK cells','Macrophages',
                                     'M1 macrophages','M2 macrophages', 
                                     'Dendritic cells','Neutrophils', 'B cells','Mast cells'))

####################################################################################################

#### OMIC DATA
AD_wide <- AD %>%
  pivot_wider(names_from = Cells_group,
              values_from = Result)

##Ajout d'un "ID" à chaque data
AD_wide <- data.frame(ID=paste0('S',seq(1,dim(AD_wide)[1])), AD_wide)

## Retrait des "coldata" pour ne garder que les data omic (=Cells_group + Result)
AD_omic_data <- AD_wide %>%
  ungroup() %>%
  dplyr::select(-Phase,-Image,-Classification) %>%
  column_to_rownames(var='ID')

AD_omic_data <- t(AD_omic_data)

##Transformation en log2
AD_omic_data <- log2(AD_omic_data+0.001)

## coldata
coldata <- AD_wide %>%
  ungroup() %>%
  dplyr::select(ID, Phase, Image, Classification)

####################################################################################################

##Design
######## multiple sample per patients 

mydesign <- model.matrix(~Classification, data = coldata)

fit <- lmFit(AD_omic_data, design = mydesign)

cm <- makeContrasts(AD_vs_T = ClassificationAD_ectopic_adenomyosis, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="AD_vs_T",n=dim(efit)[1],genelist = rownames(AD_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result


