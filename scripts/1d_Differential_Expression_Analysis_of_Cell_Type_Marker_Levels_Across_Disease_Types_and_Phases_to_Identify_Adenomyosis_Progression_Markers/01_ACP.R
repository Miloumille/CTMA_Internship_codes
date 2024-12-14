##Packages

library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)
library(openxlsx)

library(factoextra)
library(FactoMineR)

ANALYSE EN COMPOSANTE PRINCIPALE, POUR VOIR SI OUTLIER A ENLEVER
####################################################################################################

##Raw data
cells_multiplex <- read_excel('../data/cell_names_1.xlsx')

##Modified data - eutopic
eutopic_cells <- cells_multiplex %>%
  group_by(Phase,Image,Classification,Cells_subgroup) %>%
  summarise(Result=sum(Result)) %>% 
  filter(Classification %in% c("AD_eutopic", "Control_eutopic"))

##Factor --> ordre données dans graph
eutopic_cells$Classification <- factor(eutopic_cells$Classification, levels = c('Control_eutopic', 'AD_eutopic'))

#### OMIC DATA - eutopic
eutopic_wide <- eutopic_cells %>%
  pivot_wider(names_from = Cells_subgroup,
              values_from = Result)

##Ajout d'un "ID" à chaque data
eutopic_wide <- data.frame(ID=paste0('S',seq(1,dim(eutopic_wide)[1])), eutopic_wide)

## Retrait des "coldata" pour ne garder que les data omic (=Cells_subgroup + Result)
eutopic_omic_data <- eutopic_wide %>%
  ungroup() %>%
  dplyr::select(-Phase,-Image,-Classification) %>%
  column_to_rownames(var='ID')

eutopic_omic_data <- t(eutopic_omic_data)

##Transformation en log2
eutopic_omic_data <- log2(eutopic_omic_data+0.001)

## coldata
coldata <- eutopic_wide %>%
  ungroup() %>%
  dplyr::select(ID, Phase, Image, Classification)

####################################################################################################

tokeep <- colSums(is.na(eutopic_omic_data))==0

eutopic_omic_data <- eutopic_omic_data[,tokeep==T]
coldata <- coldata[tokeep==T,]

pca <- PCA(t(eutopic_omic_data),  graph = F)
fviz_pca_ind(pca, col.ind = coldata$Classification ,geom='text', repel = T ,title='PCA',palette = c('steelblue1', 'tan1'),addEllipses=TRUE,ellipse.level=0.95)


