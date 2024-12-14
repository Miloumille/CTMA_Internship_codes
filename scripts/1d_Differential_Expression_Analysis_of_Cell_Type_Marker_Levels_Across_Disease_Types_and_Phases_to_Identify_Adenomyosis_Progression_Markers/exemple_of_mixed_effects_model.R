library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyverse)
library(purrr)


###Préparation des données
data_multiplex <- read_excel('../data/CD20_CD15_CD117_PanelA_PanelB_JAMBROISE_WACHEUL.xlsx')
data_multiplex$Cells[data_multiplex$Cells=='T cells'] <- 'T_cells'


AD_ectopic_eutopic <- data_multiplex %>%
  group_by(Phase,Image,Classification,IHC,Type) %>%
  summarise(Result=sum(Result)) %>% 
  filter(Classification %in% c("AD_ectopic_adenomyosis", "AD_eutopic", "Control_eutopic")) %>%
  mutate(Classification=factor(Classification,levels=c('Control_eutopic','AD_eutopic','AD_ectopic_adenomyosis')))


table(AD_ectopic_eutopic$IHC)

#AD_ectopic_eutopic <- AD_ectopic_eutopic %>%
#  filter(is.element(IHC,c('CD117','CD15','CD163','CD20','CD3'))==T)

table(AD_ectopic_eutopic$IHC)
         
         
AD_ectopic_eutopic %>% 
  ggplot(aes(
    x = Classification,
    y = log2(Result+0.001),
    col= Classification)) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(IHC~.,scales='free')


################# 

summary_wide <- AD_ectopic_eutopic %>%
  pivot_wider(names_from = IHC,values_from = Result)
summary_wide <- data.frame(ID=paste0('S',seq(1,dim(summary_wide)[1])),summary_wide)

omic_data <- summary_wide %>%
  ungroup() %>%
  dplyr::select(-Image,-Type,-Phase,-Classification) %>%
  column_to_rownames(var='ID')
omic_data <- t(omic_data)

omic_data <- log2(omic_data+0.001)

coldata <- summary_wide %>%
  ungroup() %>%
  dplyr::select(ID,Image,Phase,Classification) %>%
  mutate(Classification=factor(Classification,levels=c('Control_eutopic','AD_eutopic','AD_ectopic_adenomyosis')))

table(coldata$Classification)

mydesign <- model.matrix(~Classification,data = coldata)

corfit <- duplicateCorrelation(omic_data,design=mydesign,block=coldata$Image)              ##
fit <- lmFit(omic_data,design=mydesign,block=coldata$Image,correlation=corfit$consensus)   ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS

cm <- makeContrasts( Adenomyosis_vs_control = ClassificationAD_ectopic_adenomyosis, levels=mydesign)   ## si on veut faire la différence entre Adenomyosis et Transition , on remplace "TypeAdenomyosis" par "TypeAdenomyosis-TypeTransition"

cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result




