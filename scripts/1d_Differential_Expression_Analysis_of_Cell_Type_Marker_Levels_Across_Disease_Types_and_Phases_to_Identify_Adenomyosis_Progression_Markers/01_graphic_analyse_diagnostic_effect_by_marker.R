library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)


data_multiplex <- read_excel('../data/CD20_CD15_CD117_PanelA_JAMBROISE_WACHEUL.xlsx')
data_multiplex$Cells[data_multiplex$Cells=='T cells'] <- 'T_cells'

data_multiplex %>%
  filter(Type == 'Adenomyosis',Cells=='B cells',Phase=='Menstrual') |> data.frame()

summary_multiplex <- data_multiplex %>%
  group_by(Image,Type,Cells,Phase,Classification) %>%
  summarise(Result=sum(Result))

summary_multiplex$Type <- factor(summary_multiplex$Type,levels=c('Control','Transition','Adenomyosis'))

pdf('../results/boxplot_diagnostic_log2.pdf',width=8,height = 12)
ggplot(summary_multiplex,mapping = aes(x = Type ,y = log2(Result+0.001),col=Type))+
  geom_boxplot()+
  geom_jitter(width = 0.1,size=0.5)+
  facet_grid(Cells~Phase,scales='free_y')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../results/boxplot_diagnostic.pdf',width=8,height = 12)
ggplot(summary_multiplex,mapping = aes(x = Type ,y = Result,col=Type))+
  geom_boxplot()+
  geom_jitter(width = 0.1,size=0.5)+
  facet_grid(Cells~Phase,scales='free_y')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


  

################## modèle avec data merged (ignore Phase variable)
  ## modèle chaque phase séparamment
   ### effet interaction entre type et phase
### refaire la même analyse sur IHC (donc sans sommer )


summary_wide <- summary_multiplex %>%
  pivot_wider(names_from = Cells,values_from = Result)
summary_wide <- data.frame(ID=paste0('S',seq(1,dim(summary_wide)[1])),summary_wide)


omic_data <- summary_wide %>%
  ungroup() %>%
  dplyr::select(-Image,-Type,-Phase,-Classification) %>%
  column_to_rownames(var='ID')
omic_data <- t(omic_data)

omic_data <- log2(omic_data+0.001)

coldata <- summary_wide %>%
  ungroup() %>%
  dplyr::select(ID,Image,Type,Phase,Classification) 


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
