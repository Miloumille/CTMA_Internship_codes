library(WriteXLS)
library(edgeR)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(WriteXLS)


counts <- read_excel('../data/Project heartmates non ischemic vs ischemic_SMIT1 implication_count tables with H49.xlsx')
colnames(counts)[-1] <- paste0('H',colnames(counts)[-1])
dge0 <- DGEList(counts=counts %>% dplyr::select(-Gene))

dge0 <- calcNormFactors(dge0, method = "TMM")
logcpm <- cpm(dge0,log = T,prior.count = 1,normalized.lib.sizes = T)
rownames(logcpm) <- counts$Gene

logcpm <- data.frame(logcpm) %>%
  rownames_to_column(var='SYMBOL')

### selection of variable genes


### selection of housekeeping genes

logcpm_selected_hskp <- logcpm %>%
  filter(is.element(SYMBOL,c('RPL32','RPLP0','B2M','GAPDH')))

###

expression_long <- logcpm_selected_hskp %>%
  pivot_longer(cols = where(is.numeric),names_to = "SampleName", values_to = "logCPM") %>%
  mutate(SampleName = fct_inorder(SampleName))

pdf('../results/Housekeeping_stability.pdf',width=10,height = 4)  
ggplot(expression_long, aes(x = SampleName, y = logCPM,col=SYMBOL)) + 
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_line(aes(group = SYMBOL), lty = 2)+
  ylim(0,15)
dev.off()




