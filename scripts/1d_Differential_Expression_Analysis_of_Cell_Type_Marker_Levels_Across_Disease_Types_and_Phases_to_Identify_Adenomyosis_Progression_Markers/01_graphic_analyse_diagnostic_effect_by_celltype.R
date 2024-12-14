library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)
library(openxlsx)

data_multiplex <- read_excel('../data/CD20_CD15_CD117_PanelA_JAMBROISE_WACHEUL.xlsx')
data_multiplex$Cells[data_multiplex$Cells=='T cells'] <- 'T_cells'

data_multiplex %>%
  filter(Type == 'Adenomyosis',Cells=='B cells',Phase=='Menstrual') |> data.frame()

summary_multiplex <- data_multiplex %>%
  group_by(Image,Type,Cells,Phase) %>%
  summarise(Result=sum(Result))

summary_multiplex$Type <- factor(summary_multiplex$Type,levels=c('Control','Transition','Adenomyosis'))

summary_multiplex <- summary_multiplex %>%
  filter(Cells!='NA',Cells!='Neutrophils')


pdf('../results/diagnostic_effect_by_celltype/scatterplot_diagnostic_phase_separated.pdf',width=6,height = 10)
ggplot(summary_multiplex,mapping = aes(x = Type ,y = Result,col=Type))+
  geom_jitter(width = 0.1,size=0.5,alpha=0.6)+
  stat_summary(geom="point",fun="mean",size=1.5,alpha=0.8)+
  facet_grid(Cells~Phase,scales='free_y')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../results/diagnostic_effect_by_celltype/scatterplot_diagnostic_phase_separated_log2.pdf',width=6,height = 10)
ggplot(summary_multiplex,mapping = aes(x = Type ,y = log2(Result+0.001),col=Type))+
  geom_jitter(width = 0.1,size=0.5,alpha=0.6)+
  stat_summary(geom="point",fun="mean",size=1.5,alpha=0.8)+
  facet_grid(Cells~Phase,scales='free_y')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf('../results/diagnostic_effect_by_celltype/scatterplot_diagnostic_phase_merged.pdf',width=4,height = 10)
ggplot(summary_multiplex,mapping = aes(x = Type ,y = Result,col=Type))+
  geom_jitter(width = 0.1,size=0.5,alpha=0.6)+
  stat_summary(geom="point",fun="mean",size=1.5,alpha=0.8)+
  facet_grid(Cells~.,scales='free_y')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../results/diagnostic_effect_by_celltype/scatterplot_diagnostic_phase_merged_log2.pdf',width=4,height = 10)
ggplot(summary_multiplex,mapping = aes(x = Type ,y = log2(Result+0.001),col=Type))+
  geom_jitter(width = 0.1,size=0.5,alpha=0.6)+
  stat_summary(geom="point",fun="mean",size=1.5,alpha=0.8)+
  facet_grid(Cells~.,scales='free_y')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



#################  MODEL WITH ALL PHASE MERGED


summary_wide <- summary_multiplex %>%
  pivot_wider(names_from = Cells,values_from = Result)
summary_wide <- data.frame(ID=paste0('S',seq(1,dim(summary_wide)[1])),summary_wide)


omic_data <- summary_wide %>%
  ungroup() %>%
  dplyr::select(-Image,-Type,-Phase,) %>%
  column_to_rownames(var='ID')
omic_data <- t(omic_data)

omic_data <- log2(omic_data+0.001)

coldata <- summary_wide %>%
  ungroup() %>%
  dplyr::select(ID,Image,Type,Phase) 


mydesign <- model.matrix(~Type,data = coldata)
fit <- lmFit(omic_data, mydesign)

cm <- makeContrasts( Adenomyosis_vs_control = TypeAdenomyosis, levels=mydesign)
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result



#################  MODEL WITH ONLY Secretory_PHASE


summary_wide <- summary_multiplex %>%
  filter(Phase=='Secretory') %>%
  pivot_wider(names_from = Cells,values_from = Result)
summary_wide <- data.frame(ID=paste0('S',seq(1,dim(summary_wide)[1])),summary_wide)


omic_data <- summary_wide %>%
  ungroup() %>%
  dplyr::select(-Image,-Type,-Phase,) %>%
  column_to_rownames(var='ID')
omic_data <- t(omic_data)

omic_data <- log2(omic_data+0.001)

coldata <- summary_wide %>%
  ungroup() %>%
  dplyr::select(ID,Image,Type,Phase) 


mydesign <- model.matrix(~Type,data = coldata)
fit <- lmFit(omic_data, mydesign)

cm <- makeContrasts( Adenomyosis_vs_control = TypeAdenomyosis, levels=mydesign)
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)


result <- topTable(efit, coef="Adenomyosis_vs_control",n=dim(efit)[1],genelist = rownames(omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

write.xlsx(result,'../results/diagnostic_effect_by_celltype/diagnostic_effect_in_secretory_phase.xlsx')





################## modèle avec data merged (ignore Phase variable)
  ## modèle chaque phase séparamment
   ### effet interaction entre type et phase
### refaire la même analyse sur IHC (donc sans sommer )



