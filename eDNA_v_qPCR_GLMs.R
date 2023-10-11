#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")

#load libraries
library(tidyverse)

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#load data (see eDNA_v_qPCR_DataProcessing.R)
AtlSalmon_data <- read.csv("qPCR_eDNA_SalmoSalar_CleanedData_10Oct2023.csv") %>% 
  select(!(X))

PinkSalmon_data <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_10Oct2023.csv") %>% 
  select(!(X))

ArcticCharr_data <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_10Oct2023.csv") %>% 
  select(!(X))

#calculate QuantMean
AtlSalmon_data <- AtlSalmon_data %>% 
  mutate(QuantMean = rowMeans(select(.,Quant1,Quant2,Quant3), na.rm=TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)

plot(log(AtlSalmon_data$QuantMean), log(AtlSalmon_data$X12Steleo))








#model variable selections

 




#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
