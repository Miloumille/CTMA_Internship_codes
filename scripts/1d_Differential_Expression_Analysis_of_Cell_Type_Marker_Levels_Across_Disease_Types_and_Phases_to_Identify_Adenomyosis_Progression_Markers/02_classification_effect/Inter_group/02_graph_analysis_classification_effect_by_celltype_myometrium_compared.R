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

##Modified data 
myometrium_cells <- cells_multiplex %>%
  group_by(Phase,Image,Classification,Cells_group) %>%
  summarise(Result=sum(Result)) %>% 
  filter(Classification %in% c("Transition_myometrium", "Control_myometrium", "AD_myometrium")) %>%
  filter(Cells_group!='NA',Cells_group!='Neutrophils')


##Factor --> ordre données dans graph
myometrium_cells$Classification <- factor(myometrium_cells$Classification, levels = c('Control_myometrium', 'Transition_myometrium', 'AD_myometrium'))

####################################################################################################

##GRAPH
#PC PORTABLE

####################################################################################################

#### OMIC DATA
myometrium_wide <- myometrium_cells %>%
  pivot_wider(names_from = Cells_group,
              values_from = Result)

##Ajout d'un "ID" à chaque data
myometrium_wide <- data.frame(ID=paste0('S',seq(1,dim(myometrium_wide)[1])), myometrium_wide)

## Retrait des "coldata" pour ne garder que les data omic (=Cells_group + Result)
myometrium_omic_data <- myometrium_wide %>%
  ungroup() %>%
  dplyr::select(-Phase,-Image,-Classification) %>%
  column_to_rownames(var='ID')

myometrium_omic_data <- t(myometrium_omic_data)

##Transformation en log2
myometrium_omic_data <- log2(myometrium_omic_data+0.001)

## coldata
coldata <- myometrium_wide %>%
  ungroup() %>%
  dplyr::select(ID, Phase, Image, Classification)

####################################################################################################

##Design 
######## unique sample per patients
mydesign <- model.matrix(~Classification, data = coldata)

fit <- lmFit(myometrium_omic_data, design = mydesign)

####################################################################################################

#################  MODEL WITH ALL PHASE MERGED Comparaison Adenomyosis_vs_control
cm <- makeContrasts(Adenomyosis_vs_control = ClassificationAD_myometrium, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(myometrium_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

###Export results
write.xlsx(result,'../results/02_classification_effect/merged_phases/cell_groups_myometrium_compared_AD_vs_control.xlsx')

####################################################################################################

#################  MODEL WITH ALL PHASE MERGED  Comparaison Transition_vs_control
cm <- makeContrasts(Transition_vs_control = ClassificationTransition_myometrium, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="Transition_vs_control",n=dim(efit)[1],genelist = rownames(myometrium_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

###Export results
write.xlsx(result,'../results/02_classification_effect/merged_phases/cell_groups_myometrium_compared_T_vs_control.xlsx')


####################################################################################################

#################  MODEL WITH ALL PHASE MERGED  Comparaison AD_vs_transition
cm <- makeContrasts(AD_vs_transition = ClassificationAD_myometrium-ClassificationTransition_myometrium, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

result <- topTable(efit, coef="AD_vs_transition",n=dim(efit)[1],genelist = rownames(myometrium_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

###Export results
write.xlsx(result,'../results/02_classification_effect/merged_phases/cell_groups_myometrium_compared_AD_vs_transition.xlsx')



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

fit <- lmFit(myometrium_omic_data, design = mydesign)

efit <- eBayes(fit,trend=T,robust = T)
efit

topTable(efit, coef="ClassificationAD_myometrium:PhaseSecretory",n=dim(efit)[1],genelist = rownames(myometrium_omic_data))
topTable(efit, coef="ClassificationAD_myometrium:PhaseProliferative",n=dim(efit)[1],genelist = rownames(myometrium_omic_data))
topTable(efit, coef="ClassificationAD_myometrium:PhaseMenstrual",n=dim(efit)[1],genelist = rownames(myometrium_omic_data))
