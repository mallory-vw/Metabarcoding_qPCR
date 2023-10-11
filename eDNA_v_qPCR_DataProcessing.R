#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_DataProcessing_workspace.RData")

#load libraries
library(tidyverse)

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")



#load all data
metadata <- read.csv("RiverSamplingMetadata_2019_2021.csv",na.strings="")

eDNA_data_raw <- read.csv("eDNA_results_PrelimCleaned_5Oct2023.csv", na.strings = "N/A") %>% 
  mutate(CorrectedDepth = as.integer(CorrectedDepth)) %>% 
  rename(RawReadsPerSample=RawDepthPerSample,
         RawReads=RawDepth,
         CorrectedReads=CorrectedDepth,
         SampleID=Sample)
  
AtlSalmon_qPCR_raw <- read.csv("qPCR_SalmoSalar_results_PrelimCleaned_5Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID)
PinkSalmon_qPCR_raw <- read.csv("qPCR_PinkSalmon_results_PrelimCleaned_5Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID)
ArcticCharr_qPCR_raw <- read.csv("qPCR_ArcticCharr_results_PrelimCleaned_10Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID)

######eDNA#####
#organize total reads per sample (reads of ALL vertebrates in the sample)
eDNA_ReadsSample <- eDNA_data_raw %>% 
  select(SampleID, Type, Marker, RawReadsPerSample) %>% 
  distinct()

#convert to wide
eDNA_ReadsSample_wide <- eDNA_ReadsSample %>% 
  pivot_wider(id_cols = c("SampleID","Type"),
              names_from = Marker,
              values_from = RawReadsPerSample)

#Field blanks, lab negatives, and unknown
eDNA_ReadsSample_Blank <- eDNA_ReadsSample %>% 
  filter(Type %in% c("Blank","Negative","Unknown"))

#Organize data to work with and calculate proportion of total reads
eDNA_data <- eDNA_data_raw %>% 
  mutate(CorrectedPercentReads = as.character(signif((CorrectedReads/RawReadsPerSample)*100, digits=4)))

#separate data by species
AtlSalmon_eDNA_data <- eDNA_data %>% 
  select(SampleID,Type,Marker,Taxon,CorrectedReads,CorrectedPercentReads) %>% 
  filter(Taxon == "Salmo salar") %>% 
  pivot_wider(id_cols = c("SampleID","Type"),
              names_from = Marker,
              values_from = c(CorrectedReads, CorrectedPercentReads)) %>% 
  relocate(CorrectedPercentReads_12Steleo, .after=CorrectedReads_12Steleo) %>% 
  relocate(CorrectedPercentReads_MIFISHU, .after=CorrectedReads_MIFISHU) %>% 
  relocate(CorrectedPercentReads_FISHE, .after=CorrectedReads_FISHE)
  
PinkSalmon_eDNA_data <- eDNA_data %>% 
  select(SampleID,Type,Marker,Taxon,CorrectedReads,CorrectedPercentReads) %>% 
  filter(Taxon == "Oncorhynchus gorbuscha") %>% 
  pivot_wider(id_cols = c("SampleID","Type"),
              names_from = Marker,
              values_from = c(CorrectedReads, CorrectedPercentReads)) %>% 
  relocate(CorrectedPercentReads_12Steleo, .after=CorrectedReads_12Steleo) %>% 
  relocate(CorrectedPercentReads_MIFISHU, .after=CorrectedReads_MIFISHU) %>% 
  relocate(CorrectedPercentReads_FISHE, .after=CorrectedReads_FISHE)

Charr_eDNA_data <- eDNA_data %>% 
  select(SampleID,Type,Marker,Taxon,CorrectedReads,CorrectedPercentReads) %>% 
  filter(Taxon == "Salvelinus alpinus") %>% 
  pivot_wider(id_cols = c("SampleID","Type"),
              names_from = Marker,
              values_from = c(CorrectedReads, CorrectedPercentReads)) %>% 
  relocate(CorrectedPercentReads_MIFISHU, .after=CorrectedReads_MIFISHU) %>% 
  relocate(CorrectedPercentReads_FISHE, .after=CorrectedReads_FISHE)


#####qPCR and eDNA#####
#Atlantic salmon#
#check and remove any rows that had inhibition
unique(AtlSalmon_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

AtlSalmon_qPCR_data <- AtlSalmon_qPCR_raw %>% 
  select(!IPCResult) %>% #remove the IPCResults column
  filter(Result %in% c("Detected","Suspected")) #keep only the Detected and Suspected samples

#join data, keeping ONLY the samples in both qPCR and eDNA results
AtlSalmon_data <- AtlSalmon_eDNA_data %>% 
  inner_join(AtlSalmon_qPCR_data) %>% 
  left_join(metadata)

#export
write.csv(AtlSalmon_data, "qPCR_eDNA_SalmoSalar_CleanedData_11Oct2023.csv")

#Pink salmon#
#check and remove any rows that had inhibition
unique(PinkSalmon_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

PinkSalmon_qPCR_data <- PinkSalmon_qPCR_raw %>% 
  select(!IPCResult) %>% #remove the IPCResults column
  filter(Result %in% c("Detected","Suspected")) #keep only the Detected and Suspected samples

#join data, keeping ONLY the samples in both qPCR and eDNA results
PinkSalmon_data <- PinkSalmon_eDNA_data %>% 
  inner_join(PinkSalmon_qPCR_data) %>% 
  left_join(metadata)

#export
write.csv(PinkSalmon_data, "qPCR_eDNA_PinkSalmon_CleanedData_11Oct2023.csv")

#Arctic charr#
#check and remove any rows that had inhibition
unique(ArcticCharr_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

ArcticCharr_qPCR_data <- ArcticCharr_qPCR_raw %>% 
  select(!IPCResult) %>% #remove the IPCResults column
  filter(Result %in% c("Detected","Suspected")) %>% #keep only the Detected and Suspected samples
  filter(!SampleID %in% c("GLPR","SWA2019C2","MBB2019C2")) %>%  #remove problem samples based on the Notes field
  mutate(Ct2 = replace(Ct2, which(SampleID %in% c("TOM2019C3","PAN2019C1","PAN2019C3","SWA2019C3")), NA)) %>% #replace empty replicates with NA based on Notes field
  mutate(Ct3 = replace(Ct3, which(SampleID %in% c("TOM2019C3","PAN2019C1","PAN2019C3","SWA2019C3","SWA2019C1","KIN2019C1")), NA)) #replace empty replicates with NA based on Notes field

#join data, keeping ONLY the samples in both qPCR and eDNA results
ArcticCharr_data <- Charr_eDNA_data %>% 
  inner_join(ArcticCharr_qPCR_data) %>% 
  left_join(metadata) %>% 
  filter()

#export
write.csv(ArcticCharr_data, "qPCR_eDNA_ArcticCharr_CleanedData_11Oct2023.csv")


#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_DataProcessing_workspace.RData")
