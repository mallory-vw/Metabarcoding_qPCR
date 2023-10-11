#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")

#load libraries
library(tidyverse)

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#load data (see eDNA_v_qPCR_DataProcessing.R)
AtlSalmon_data <- read.csv("qPCR_eDNA_SalmoSalar_CleanedData_11Oct2023.csv") %>% 
  select(!(X)) %>% 
  mutate(CorrectedPercentReads_12Steleo = as.character(signif(CorrectedPercentReads_12Steleo, digits=4))) %>% 
  mutate(CorrectedPercentReads_MIFISHU = as.character(signif(CorrectedPercentReads_MIFISHU, digits=4))) %>% 
  mutate(CorrectedPercentReads_FISHE = as.character(signif(CorrectedPercentReads_FISHE, digits=4)))

PinkSalmon_data <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_11Oct2023.csv") %>% 
  select(!(X)) %>% 
  mutate(CorrectedPercentReads_12Steleo = as.character(signif(CorrectedPercentReads_12Steleo, digits=4))) %>% 
  mutate(CorrectedPercentReads_MIFISHU = as.character(signif(CorrectedPercentReads_MIFISHU, digits=4))) %>% 
  mutate(CorrectedPercentReads_FISHE = as.character(signif(CorrectedPercentReads_FISHE, digits=4)))

ArcticCharr_data <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_11Oct2023.csv") %>% 
  select(!(X)) %>% 
  mutate(CorrectedPercentReads_MIFISHU = as.character(signif(CorrectedPercentReads_MIFISHU, digits=4))) %>% 
  mutate(CorrectedPercentReads_FISHE = as.character(signif(CorrectedPercentReads_FISHE, digits=4)))

#calculate QuantMean
AtlSalmon_data <- AtlSalmon_data %>% 
  mutate(QuantMean = rowMeans(select(.,Quant1,Quant2,Quant3), na.rm=TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)

#normalize CorrectedPercentReads by Vol sampled
#normalize QuantMean by DNA concentration? or also by vol filtered, since DNA concentration will also be related to vol filtered
#or just control for volume filtered in the model?

test <- AtlSalmon_data %>% 
  select(SampleID,CorrectedPercentReads_12Steleo, VolFiltered, QuantMean, DNAConc_pg_uL) %>% 
  mutate_at(vars(CorrectedPercentReads_12Steleo,DNAConc_pg_uL),as.numeric) %>% 
  mutate(NormPercReads = CorrectedPercentReads_12Steleo/VolFiltered,
         NormQuantMean = QuantMean/VolFiltered,
         DNANormQuantMean = QuantMean/DNAConc_pg_uL)

testlm <- glm(test$NormPercReads ~ test$NormQuantMean)

anova(testlm)
summary(testlm)

#model variable selections

 




#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
