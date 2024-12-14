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
cells_multiplex <- read_excel('../data/cell_names_2.xlsx')

##Modified data - eutopic
eutopic_cells <- cells_multiplex %>%
  group_by(Phase,Image,Classification,Cells_group) %>%
  summarise(Result=sum(Result)) %>% 
  filter(Classification %in% c("AD_eutopic", "Control_eutopic", "Transition_eutopic"))

eutopic_cells <- eutopic_cells %>%
  filter(Cells_group!='NA',Cells_group!='Neutrophils')

##Factor --> ordre données dans graph
eutopic_cells$Classification <- factor(eutopic_cells$Classification, levels = c('Control_eutopic','Transition_eutopic', 'AD_eutopic'))

####################################################################################################

##GRAPH
##SCRIPT PC PORTABLE

AD$Cells_subgroup <- factor(AD$Cells_subgroup,
                            levels=c('T cells',
                                     'Cytotoxic T cells','Helper 1 T cells',
                                     'Helper 2 T cells','NK cells','Macrophages',
                                     'M1 macrophages','M2 macrophages', 
                                     'Dendritic cells','Neutrophils', 'B cells','Mast cells'))
####################################################################################################

#### OMIC DATA - eutopic
eutopic_wide <- eutopic_cells %>%
  pivot_wider(names_from = Cells_group,
              values_from = Result)

##Ajout d'un "ID" à chaque data
eutopic_wide <- data.frame(ID=paste0('S',seq(1,dim(eutopic_wide)[1])), eutopic_wide)

## Retrait des "coldata" pour ne garder que les data omic (=Cells_group + Result)
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

##Design - eutopic
######## unique sample per patients
mydesign <- model.matrix(~Classification, data = coldata)

fit <- lmFit(eutopic_omic_data, design = mydesign)

####################################################################################################

#################  MODEL WITH ALL PHASE MERGED Comparaison Adenomyosis_vs_control
cm <- makeContrasts(Adenomyosis_vs_control = ClassificationAD_eutopic, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(eutopic_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/02_classification_effect/merged_phases/cell_groups_eutopic_AD_versus_control.xlsx')


####################################################################################################

#################  MODEL WITH ALL PHASE MERGED  Comparaison Transition_vs_control
cm <- makeContrasts(Transition_vs_control = ClassificationTransition_eutopic, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="Transition_vs_control",n=dim(efit)[1],genelist = rownames(eutopic_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/02_classification_effect/merged_phases/cell_groups_eutopic_T_versus_control.xlsx')


####################################################################################################

#################  MODEL WITH ALL PHASE MERGED  Comparaison AD_vs_transition
cm <- makeContrasts(AD_vs_transition = ClassificationAD_eutopic-ClassificationTransition_eutopic, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="AD_vs_transition",n=dim(efit)[1],genelist = rownames(eutopic_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/02_classification_effect/merged_phases/cell_groups_eutopic_AD_versus_T.xlsx')


####################################################################################################
####################################################################################################
### effet interaction entre type et phase

##Design
##choisir quel est le niveau de "référence" auquel seront comparés les autres niveaux 
#Si secretory
coldata$Phase <- factor(coldata$Phase,levels=c('Secretory', 'Proliferative','Menstrual'))
#Si Proliferative
coldata$Phase <- factor(coldata$Phase,levels=c('Proliferative','Secretory','Menstrual'))
#Si Menstrual
coldata$Phase <- factor(coldata$Phase,levels=c('Menstrual', 'Secretory','Proliferative'))


mydesign <- model.matrix(~Classification+Phase+Phase*Classification, data = coldata)

fit <- lmFit(eutopic_omic_data, design = mydesign)

efit <- eBayes(fit,trend=T,robust = T)
efit

topTable(efit, coef="ClassificationAD_eutopic:PhaseSecretory",n=dim(efit)[1],genelist = rownames(eutopic_omic_data))
topTable(efit, coef="ClassificationAD_eutopic:PhaseProliferative",n=dim(efit)[1],genelist = rownames(eutopic_omic_data))
topTable(efit, coef="ClassificationAD_eutopic:PhaseMenstrual",n=dim(efit)[1],genelist = rownames(eutopic_omic_data))
