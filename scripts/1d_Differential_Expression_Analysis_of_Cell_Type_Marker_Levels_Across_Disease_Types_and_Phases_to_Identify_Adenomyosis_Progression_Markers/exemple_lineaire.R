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
cells_multiplex <- read_excel('../data/cell_names_1_AT_grouped.xlsx')

##Modified data
AD <- cells_multiplex %>%
  group_by(Phase,Image,Classification,Cells_subgroup) %>%
  summarise(Result=sum(Result)) %>%
  filter(is.element(Classification,c('AD_eutopic','AD_ectopic_adenomyosis','AD_ectopic_transition', 'AD_ectopic_transition_dilated')))

Classification_num <- numeric()
Classification_num[AD$Classification=='AD_eutopic'] <- 0
Classification_num[AD$Classification=='AD_ectopic_transition'] <- 2
Classification_num[AD$Classification=='AD_ectopic_transition_dilated'] <- 2
Classification_num[AD$Classification=='AD_ectopic_adenomyosis'] <- 2.5

AD <- data.frame(AD,Classification_num)


pdf('../results/02_classification_effect/linear_regression.pdf', width=15,height = 25)
ggplot(AD, mapping = aes(
  x = Classification_num,
  y = log2(Result+0.001))) +
  geom_jitter(width = 0.1, size=0.5, alpha=0.6) +
  facet_wrap(~Cells_subgroup,scales='free_y',ncol=3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_smooth(method=lm)
dev.off()



#### OMIC DATA - AD

AD_wide <- AD %>%
  pivot_wider(names_from = Cells_subgroup,
              values_from = Result)

##Ajout d'un "ID" à chaque data
AD_wide <- data.frame(ID=paste0('S',seq(1,dim(AD_wide)[1])), AD_wide)

## Retrait des "coldata" pour ne garder que les data omic (=Cells_subgroup + Result)
AD_omic_data <- AD_wide %>%
  ungroup() %>%
  dplyr::select(-Phase,-Image,-Classification,-Classification_num) %>%
  column_to_rownames(var='ID')

AD_omic_data <- t(AD_omic_data)

##Transformation en log2
AD_omic_data <- log2(AD_omic_data+0.001)

## coldata
coldata <- AD_wide %>%
  ungroup() %>%
  dplyr::select(ID, Phase, Image, Classification,Classification_num)

####################################################################################################

mydesign <- model.matrix(~Classification_num, data = coldata)
corfit <- duplicateCorrelation(AD_omic_data, design = mydesign, block = coldata$Image)  ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS
fit <- lmFit(AD_omic_data, design = mydesign, block = coldata$Image, correlation = corfit$consensus)   ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS
efit <- eBayes(fit,trend=T,robust = T)
efit

topTable(efit, coef="Classification_num",n=dim(efit)[1],genelist = rownames(AD_omic_data))



