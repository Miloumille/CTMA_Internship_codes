library(ggplot2)
library(readxl)
library(dplyr)
library(limma)
library(tidyr)
library(purrr)


data_multiplex <- read_excel('../data/CD20_CD15_CD117_PanelA_JAMBROISE_WACHEUL.xlsx')
data_multiplex$Cells[data_multiplex$Cells=='T cells'] <- 'T_cells'


data_multiplex %>%
  filter(Cells=='T_cells',Phase=='Menstrual') %>%
  arrange(Classification) |> as.data.frame() 

data_summary <- data_multiplex %>%
  group_by(Phase,Image,Classification,Cells) %>%
  summarise(Result=sum(Result))

dim(data_multiplex)
dim(data_summary)


data_summary  <- data.frame(data_summary,diagnostic=map_chr(.x = data_summary$Classification,.f = function(x) strsplit(x,split='_')[[1]][1]))



pdf('../results/boxplot_classification.pdf',width=15,height = 15)
ggplot(data_summary,mapping = aes(x = Classification ,y = Result,col=Classification))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(Cells~diagnostic,scales='free')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


  


summary_multiplex <- data_multiplex %>%
  group_by(Image,Type,Cells,Phase,Classification) %>%
  summarise(Result=sum(Result))


