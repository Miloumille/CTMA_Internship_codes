library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)
library(openxlsx)


##Raw data
cells_multiplex <- read_excel('../data/cell_names_2.xlsx')

summary_multiplex_cells <- cells_multiplex %>%
  group_by(Image,Type,Cells_group,Phase) %>%
  summarise(Result=sum(Result))

summary_multiplex_cells$Type <- factor(summary_multiplex_cells$Type,levels=c('Control','Transition','Adenomyosis'))

summary_multiplex_cells <- summary_multiplex_cells %>%
  filter(Cells_group!='NA',Cells_group!='Neutrophils')



##GRAPH
pdf('../results/01_diagnostic_effect/scatterplot_diagnostic_cells_group_phase_separated_log2.pdf',width=6,height = 10)
ggplot(summary_multiplex_cells,mapping = aes(
  x = Type,
  y = log2(Result+0.001),
  col=Type)) +
  geom_jitter(width = 0.1, size=0.5, alpha=0.6) +
  stat_summary(geom="point",fun="mean",size=1.5,alpha=0.8) +
  facet_grid(Cells_group~Phase,scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf('../results/01_diagnostic_effect/scatterplot_diagnostic_cells_group_phase_merged_log2.pdf',width=4,height = 10)
ggplot(summary_multiplex_cells,mapping = aes(
  x = Type,
  y = log2(Result+0.001),
  col=Type)) +
  geom_jitter(width = 0.1,size=0.5,alpha=0.6) +
  stat_summary(geom="point",fun="mean",size=1.5,alpha=0.8)+
  facet_grid(Cells_group~.,scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


##Coldata
summary_wide <- summary_multiplex_cells %>%
  pivot_wider(names_from = Cells_group,values_from = Result)
summary_wide <- data.frame(ID=paste0('S',seq(1,dim(summary_wide)[1])),summary_wide)

omic_data <- summary_wide %>%
  ungroup() %>%
  dplyr::select(-Image,-Type,-Phase) %>%
  column_to_rownames(var='ID')
omic_data <- t(omic_data)

omic_data <- log2(omic_data+0.001)

coldata <- summary_wide %>%
  ungroup() %>%
  dplyr::select(ID,Image,Type,Phase) 


mydesign <- model.matrix(~Type,data = coldata)
fit <- lmFit(omic_data, mydesign)


#################  MODEL WITH ALL PHASE MERGED - AD versus control
cm <- makeContrasts(Adenomyosis_vs_control = TypeAdenomyosis, levels=mydesign)   ## si on veut faire la différence entre Adenomyosis et Transition , on remplace "TypeAdenomyosis" par "TypeAdenomyosis-TypeTransition"
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/01_diagnostic_effect/diagnostic_cells_group_phase_merged_AD_versus_control.xlsx')



#################  MODEL WITH ALL PHASE MERGED - Transition versus control
cm <- makeContrasts(Transition_vs_control = TypeTransition, levels=mydesign)   ## si on veut faire la différence entre Adenomyosis et Transition , on remplace "TypeAdenomyosis" par "TypeAdenomyosis-TypeTransition"
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


result <- topTable(efit, coef="Transition_vs_control",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/01_diagnostic_effect/diagnostic_cells_group_phase_merged_T_versus_control.xlsx')



#################  MODEL WITH ALL PHASE MERGED - AD versus Transition
cm <- makeContrasts(Adenomyosis_versus_transition = TypeAdenomyosis-TypeTransition, levels=mydesign)   ## si on veut faire la différence entre Adenomyosis et Transition , on remplace "TypeAdenomyosis" par "TypeAdenomyosis-TypeTransition"
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


result <- topTable(efit, coef="Adenomyosis_versus_transition",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/01_diagnostic_effect/diagnostic_cells_group_phase_merged_AD_versus_T.xlsx')



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


mydesign <- model.matrix(~Type+Phase+Phase*Type, data = coldata)

fit <- lmFit(omic_data, design = mydesign)

efit <- eBayes(fit,trend=T,robust = T)
efit

topTable(efit, coef="TypeAdenomyosis:PhaseSecretory",n=dim(efit)[1],genelist = rownames(omic_data))
topTable(efit, coef="TypeAdenomyosis:PhaseProliferative",n=dim(efit)[1],genelist = rownames(omic_data))
topTable(efit, coef="TypeAdenomyosis:PhaseMenstrual",n=dim(efit)[1],genelist = rownames(omic_data))
