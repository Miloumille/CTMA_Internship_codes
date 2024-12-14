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
  group_by(Phase,Image,Classification_grouped,IHC) %>%
  summarise(Result=sum(Result)) %>%
  filter(Classification_grouped %in% c("AD_ectopic", "AD_eutopic"))


##Factor --> ordre données dans graph
AD$Classification_grouped <- factor(AD$Classification_grouped, levels = c('AD_eutopic', 'AD_ectopic'))

####################################################################################################

##GRAPH
pdf('../results/02_classification_effect/merged_phases/AD_vs_eutopic_IHC.pdf', width=15,height = 25)
ggplot(AD, mapping = aes(
  x = Classification_grouped,
  y = log2(Result+0.001),
  col = Classification_grouped)) +
  geom_jitter(width = 0.1, size=0.5, alpha=0.6) +
  facet_grid(IHC~.,scales='free_y') + 
  stat_summary(geom="point", fun="mean", size=1.5, alpha=0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


####################################################################################################

#### OMIC DATA
AD_wide <- AD %>%
  pivot_wider(names_from = IHC,
              values_from = Result)

##Ajout d'un "ID" à chaque data
AD_wide <- data.frame(ID=paste0('S',seq(1,dim(AD_wide)[1])), AD_wide)

## Retrait des "coldata" pour ne garder que les data omic (=IHC + Result)
AD_omic_data <- AD_wide %>%
  ungroup() %>%
  dplyr::select(-Phase,-Image,-Classification_grouped) %>%
  column_to_rownames(var='ID')

AD_omic_data <- t(AD_omic_data)

##Transformation en log2
AD_omic_data <- log2(AD_omic_data+0.001)

## coldata
coldata <- AD_wide %>%
  ungroup() %>%
  dplyr::select(ID, Phase, Image, Classification_grouped)

####################################################################################################

##Design
######## multiple sample per patients [Ici : Classification_grouped]

mydesign <- model.matrix(~Classification_grouped, data = coldata)

corfit <- duplicateCorrelation(AD_omic_data, design = mydesign, block = coldata$Image)  ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS
fit <- lmFit(AD_omic_data, design = mydesign, block = coldata$Image, correlation = corfit$consensus)   ##  SPEFICIF FOR ANALYSIS OF DATASET WITH MULTIPLE SAMPLES PER PATIENTS

####################################################################################################

##Comparaison AD_vs_eutopic
cm <- makeContrasts(AD_vs_eutopic = Classification_groupedAD_ectopic, levels=mydesign)   
cm
fit <- contrasts.fit(fit, cm)
efit <- eBayes(fit,trend=T,robust = T)

####################################################################################################

####Resultats avec pvalue
result <- topTable(efit, coef="AD_vs_eutopic",n=dim(efit)[1],genelist = rownames(AD_omic_data))
result <- result %>% select(ID,logFC,P.Value,adj.P.Val) %>%
  mutate(logFC=round(logFC,2))

result

####################################################################################################

###Export results
write.xlsx(result,'../results/02_classification_effect/merged_phases/AD_vs_eutopic_IHC.xlsx')

####################################################################################################

##Volcanoplot

# add a column of NAs
result$diffexpressed <- "NO"

result
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result$diffexpressed[result$logFC > 0.6 & result$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result$diffexpressed[result$logFC < -0.6 & result$P.Value < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" 
result$delabel <- NA
result$delabel[result$diffexpressed != "NO"] <- result$ID[result$diffexpressed != "NO"]


p <- ggplot(result, mapping = aes(
  x = logFC,
  y = -log10(P.Value),
  col=diffexpressed, 
  label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


pdf('../results/02_classification_effect/merged_phases/volcanoplot.pdf')
p2
dev.off()
