#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabarcoding_qPCR/eDNA_v_qPCR_DataProcessing_workspace.RData")

#load libraries
library("tidyverse")
library("ggmap")
library("metR")
library("tmap")
library("tmaptools")
library("RColorBrewer")
library("janitor")

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load all raw data#####
sample_metadata_raw <- read.csv("Metadata/SampleMetadata_EnviroData_7Nov2023.csv",na.strings=c("","NA")) %>% #NOTE: Date edited in Excel to Month/Day/Year
  rename(Name=RiverName,
         Code=RiverCode,
         Verified=Sample.Verified) %>% 
  mutate(Year=factor(Year),
         VolFil=as.numeric(VolFil)) 
  
site_metadata_raw <- read.csv("Metadata/SiteMetadata_EnviroData_7Nov2023.csv",na.strings=c("","NA")) %>% #NOTE: Date as Month/Day/Year
  rename(Name=RiverName,
         Code=RiverCode) %>% 
  mutate(Year=factor(Year)) 


eDNA_data_raw <- read.csv("eDNAData/eDNA_results_PrelimCleaned_8Nov2023.csv", na.strings = c("N/A","NA","")) %>% 
  dplyr::select(!X) %>% 
  mutate(CorrectedVertReadsPerSample=as.integer(CorrectedVertReadsPerSample),
         CorrectedReads=as.integer(CorrectedReads),
         Type=factor(Type))


AtlSalmon_qPCR_raw <- read.csv("qPCRData/qPCR_AtlSalmon_results_PrelimCleaned_5Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID,
         DNAConc=DNAConc_pg_uL)
PinkSalmon_qPCR_raw <- read.csv("qPCRData/qPCR_PinkSalmon_results_PrelimCleaned_5Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID,
         DNAConc=DNAConc_pg_uL)
ArcticCharr_qPCR_raw <- read.csv("qPCRData/qPCR_ArcticCharr_results_PrelimCleaned_10Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID,
         DNAConc=DNAConc_pg_uL)

#####Organize Sample Site data#####

#Fix sample names and site names mismatch - connect by Code
Samples <- sample_metadata_raw %>%  #Pull out unique Sample Code/Name combo (176 codes)
  dplyr::select(Name,Code) %>% 
  rename(Sample=Name) %>% 
  unique()

Sites <- site_metadata_raw %>%  #Pull out unique Site Code/Name combo (178 codes)
  dplyr::select(Name,Code) %>% 
  rename(Site=Name) %>% 
  unique()

Names <- Sites %>%  #Join sample and Sites
  full_join(Samples, by="Code")

Names_Issues <- Names %>%  #Figure out names with Issues
  group_by(Code) %>% 
  mutate(Count=n(), #Determine duplicated codes (different names)
         Match=identical(Site,Sample)) %>% #Determine mismatches between Site and Sample Name (same Code)
  filter(Count > 1 | Match == FALSE) 

Names <- Names %>% 
  mutate(NewName=case_when(Code == "FOR" ~ "Forteau Brook",#Manually correct the codes with more than one name/spelling
                           Code == "IKA" ~ "Ikarut River",
                           Code == "KIN" ~ "Kingurutik River",
                           Code == "PAN" ~ "Pangertok River",
                           Code == "PIN" ~ "Pinware River",
                           Code == "SAN" ~ "Sandhill River",
                           Code == "CHR" ~ "St. Charles River",
                           Code == "SUS" ~ "Susan River",
                           Code == "HUR2" ~ "Hunt River 20 km upstream",
                           Code == "HUR3" ~ "Hunt River furthest upstream",
                           Code == "IKI" ~ "Ikinet",
                           Code == "MLR" ~ "Mulligan River",
                           Code == "NAS" ~ "Naskapi River",
                           Code == "PAM" ~ "Pamialuik",
                           Code == "PIN2" ~ "Pinware miidway upriver",
                           Code == "PIN3" ~ "Pinware top watershed site",
                           Code == "PRB" ~ "Southwest Brook Paradise River",
                           Code == "TOM2" ~ "Tom Luscomb middle site",
                           Code == "TOM3" ~ "Tom Luscomb top site",
                           Code == "CMP" ~ "Campbelton River",
                           Code == "CAR" ~ "Cat Arms River",
                           Code == "FRA" ~ "Fraser River",
                           Code == "GOS" ~ "Goose River",
                           Code == "MIC" ~ "Michaels River",
                           Code == "PPB" ~ "Partridge Pond Brook",
                           Code == "TSP" ~ "Traverspine River",
                           TRUE ~ NA),
         CorrectName = coalesce(NewName,Site)) %>% #combine the corrected names with the original name list
  dplyr::select(!c(Site,Sample,NewName)) %>% #remove extraneous columns
  unique() #remove duplicates in Code and Name combo

#join all Sample data and Site data with corrected names
sample_metadata <- sample_metadata_raw %>% 
  left_join(Names) %>% 
  left_join(site_metadata_raw, by=c("Code","Year","Date")) %>%  #merge in site data, manually changed Crooked River, Shinney's River side channel codes to match properly
  rename(SampleName = Name.x, #rename fields
         SiteName = Name.y,
         Name=CorrectName,
         SampleLat = Latitude.x,
         SampleLon = Longitude.x,
         SiteLat = Latitude.y,
         SiteLon = Longitude.y,
         SampleNotes = Notes.x,
         SiteNotes = Notes.y,
         FullTrans=FullTransect.y.n,
         MeanDepth=Mean.Depth,
         MeanFlow=Mean.Flow,
         StartTime=StartTime.NDT,
         EndTime=EndTime.NDT) %>% 
  separate(col=Date, into=c("Month","Day","Year2"),sep="/", remove=F) %>% #Split Date into separate Month/Day columns
  mutate(Month=case_when(Month==8 ~ "August", #Replace Month number with name and set factor levels
                         Month==9 ~ "September",
                         Month==10 ~ "October",
                         Month==11 ~ "November",
                         Month==12 ~ "December"),
         Month=factor(Month, levels=c("August","September","October","November","December")),
         Date=mdy(Date)) %>% #convert date to Date format
  mutate(NewLatPrep=case_when(Type == "Center" ~ Center.Latitude, #consolidate latitude from site to sample based on sample type, either Center, Left, or Right
                              Type == "Right" ~ Right.Latitude,
                              Type == "Left" ~ Left.Latitude,
                              Type == "Negative" ~ Center.Latitude),
         NewLat=case_when(is.na(NewLatPrep) & is.na(SiteLat) & is.na(SampleLat) & is.na(Left.Latitude) & is.na(Center.Latitude) ~ Right.Latitude,
                          is.na(NewLatPrep) & is.na(SiteLat) & is.na(SampleLat) & is.na(Left.Latitude) ~ Center.Latitude,
                          is.na(NewLatPrep) & is.na(SiteLat) & is.na(SampleLat) ~ Left.Latitude,
                          is.na(NewLatPrep) & is.na(SiteLat) ~ SampleLat,
                          is.na(NewLatPrep) & is.na(SiteLat) ~ SampleLat,
                          is.na(NewLatPrep) ~ SiteLat,
                          TRUE ~ NewLatPrep)) %>% 
  mutate(NewLonPrep=case_when(Type == "Center" ~ Center.Longitude, #consolidate longitude from site to sample based on sample type, either Center, Left, or Right
                              Type == "Right" ~ Right.Longitude,
                              Type == "Left" ~ Left.Longitude,
                              Type == "Negative" ~ Center.Longitude),
       NewLon=case_when(is.na(NewLonPrep) & is.na(SiteLon) & is.na(SampleLon) & is.na(Left.Longitude) & is.na(Center.Longitude) ~ Right.Longitude,
                        is.na(NewLonPrep) & is.na(SiteLon) & is.na(SampleLon) & is.na(Left.Longitude) ~ Center.Longitude,
                        is.na(NewLonPrep) & is.na(SiteLon) & is.na(SampleLon) ~ Left.Longitude,
                        is.na(NewLonPrep) & is.na(SiteLon) ~ SampleLon,
                        is.na(NewLonPrep) & is.na(SiteLon) ~ SampleLon,
                        is.na(NewLonPrep) ~ SiteLon,
                        TRUE ~ NewLonPrep)) %>% 
  dplyr::select(!c(SampleName,SiteName,Year2,NewLatPrep,NewLonPrep,SampleLat,SampleLon,SiteLat,SiteLon,
                   Left.Latitude,Left.Longitude,Center.Latitude,Center.Longitude,Right.Latitude,Right.Longitude)) %>%  #remove extraneous columns
  relocate(c(Month,Day),.after=Year) %>% 
  relocate(c(NewLat,NewLon), .before=WPT) %>% 
  relocate(c(MeanDepth,MeanFlow), .after=River.Width) %>% 
  relocate(c(SampleNotes,SiteNotes), .before=Depth1) %>% 
  relocate(FullTrans, .after=WC) %>% 
  relocate(Name, .before=Code)

rm(Samples,Sites,Names,Names_Issues)

#get some data summaries
sample_metadata %>% 
  group_by(Year) %>% 
  summarize(Rivers = length(unique(Code)),
            Samples = length(unique(SampleID)))
            


######eDNA data#####
#organize total reads per sample (reads of ALL vertebrates in the sample)
eDNA_ReadsSample <- eDNA_data_raw %>% 
  dplyr::select(SampleID, Type, Marker, RawVertReadsPerSample) %>% 
  distinct()

#convert to wide
eDNA_ReadsSample_wide <- eDNA_ReadsSample %>% 
  pivot_wider(id_cols = c("SampleID","Type"),
              names_from = Marker,
              values_from = RawVertReadsPerSample)

#Field blanks, lab negatives, and unknown
eDNA_ReadsSample_Blank <- eDNA_ReadsSample %>% 
  filter(Type %in% c("Blank","Negative","Unknown"))

#Organize data to work with and calculate proportion of total reads, removing blanks and negatives
eDNA_data <- eDNA_data_raw %>% 
  filter(!Type %in% c("Blank","Negative")) %>%  #remove blanks and negatives
  filter(!(Type == "Unknown" & is.na(RawVertReadsPerSample))) %>% #remove the "Unknowns" that also have NA in the Total Vert Reads
  filter(RawVertReadsPerSample > 0) %>% #remove anything that had no vert reads at all
  mutate(CorrectedPropReads = as.character(signif((CorrectedReads/CorrectedVertReadsPerSample), digits=4))) %>% 
  mutate(RawPropReads = as.character(signif((RawReads/RawVertReadsPerSample), digits=4))) %>% 
  relocate(RawPropReads, .after=RawReads)

#separate data by species
AtlSalmon_eDNA_data <- eDNA_data %>% 
  filter(Taxon == "Salmo salar")

PinkSalmon_eDNA_data <- eDNA_data %>% 
  filter(Taxon == "Oncorhynchus gorbuscha")

Charr_eDNA_data <- eDNA_data %>% 
  filter(Taxon == "Salvelinus alpinus")


#####combine sample site, qPCR, and eDNA#####
#Atlantic salmon#
#check and remove any rows that had inhibition
unique(AtlSalmon_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

AtlSalmon_qPCR_data <- AtlSalmon_qPCR_raw %>% 
  dplyr::select(!IPCResult) %>% #remove the IPCResults column
  mutate(QuantMean=rowMeans(dplyr::select(., c("Quant1","Quant2","Quant3")), na.rm = TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)

#data summary
AtlSalmon_qPCR_data %>% 
  group_by(Sample_Year) %>% 
  summarize(qPCRSamples = length(unique(SampleID)))

#join data, keeping ONLY the samples in both qPCR and eDNA results
AtlSalmon_data <- AtlSalmon_eDNA_data %>% 
  inner_join(AtlSalmon_qPCR_data) %>% 
  left_join(sample_metadata) %>% 
  relocate(VolFil, .before=Sample_Year) %>% 
  relocate(RawPropReads, .after=RawReads) %>% 
  relocate(CorrectedPropReads, .after = CorrectedReads) %>% 
  relocate(c(Name,Code,ReplicateSet,Year,Month,Day,Date,NewLat,NewLon),.after=Type) %>% 
  select(!Sample_Year)
#data summary
AtlSalmon_data %>% 
  group_by(Year) %>% 
  summarize(qPCRSamples = length(unique(SampleID)))

#export
write.csv(AtlSalmon_data, "qPCR_eDNA_AtlSalmon_CleanedData_7Nov2023.csv")

#Pink salmon#
#check and remove any rows that had inhibition
unique(PinkSalmon_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

PinkSalmon_qPCR_data <- PinkSalmon_qPCR_raw %>% 
  dplyr::select(!IPCResult) %>% #remove the IPCResults column
  mutate(QuantMean=rowMeans(dplyr::select(., c("Quant1","Quant2","Quant3")), na.rm = TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)

#data summary
PinkSalmon_qPCR_data %>% 
  group_by(Sample_Year) %>% 
  summarize(qPCRSamples = length(unique(SampleID)))

#join data, keeping ONLY the samples in both qPCR and eDNA results
PinkSalmon_data <- PinkSalmon_eDNA_data %>% 
  inner_join(PinkSalmon_qPCR_data) %>% 
  left_join(sample_metadata) %>% 
  relocate(VolFil, .before=Sample_Year) %>% 
  relocate(RawPropReads, .after=RawReads) %>% 
  relocate(CorrectedPropReads, .after = CorrectedReads) %>% 
  relocate(c(Name,Code,ReplicateSet,Year,Month,Day,Date,NewLat,NewLon),.after=Type) %>% 
  dplyr::select(!Sample_Year)
#data summary
PinkSalmon_data %>% 
  group_by(Year) %>% 
  summarize(qPCRSamples = length(unique(SampleID)))
#export
write.csv(PinkSalmon_data, "qPCR_eDNA_PinkSalmon_CleanedData_7Nov2023.csv")

#Arctic charr#
#check and remove any rows that had inhibition
unique(ArcticCharr_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

ArcticCharr_qPCR_data <- ArcticCharr_qPCR_raw %>% 
  dplyr::select(!IPCResult) %>% #remove the IPCResults column
  mutate(QuantMean=rowMeans(dplyr::select(., c("Quant1","Quant2","Quant3")), na.rm = TRUE)) %>% 
  relocate(QuantMean, .after=Quant3) %>% 
  filter(!SampleID %in% c("GLPR","SWA2019C2","MBB2019C2")) %>%  #remove problem samples based on the Notes field
  mutate(Ct2 = replace(Ct2, which(SampleID %in% c("TOM2019C3","PAN2019C1","PAN2019C3","SWA2019C3")), NA)) %>% #replace empty replicates with NA based on Notes field
  mutate(Ct3 = replace(Ct3, which(SampleID %in% c("TOM2019C3","PAN2019C1","PAN2019C3","SWA2019C3","SWA2019C1","KIN2019C1")), NA)) #replace empty replicates with NA based on Notes field

#data summary
ArcticCharr_qPCR_data %>% 
  group_by(Sample_Year) %>% 
  summarize(qPCRSamples = length(unique(SampleID)))

#join data, keeping ONLY the samples in both qPCR and eDNA results
ArcticCharr_data <- Charr_eDNA_data %>% 
  inner_join(ArcticCharr_qPCR_data) %>% 
  left_join(sample_metadata) %>% 
  relocate(VolFil, .before=Sample_Year) %>% 
  relocate(RawPropReads, .after=RawReads) %>% 
  relocate(CorrectedPropReads, .after = CorrectedReads) %>% 
  relocate(c(Name,Code,ReplicateSet,Year,Month,Day,Date,NewLat,NewLon),.after=Type) %>% 
  dplyr::select(!Sample_Year)

#data summary
ArcticCharr_data %>% 
  group_by(Year) %>% 
  summarize(qPCRSamples = length(unique(SampleID)))

#export
write.csv(ArcticCharr_data, "qPCR_eDNA_ArcticCharr_CleanedData_7Nov2023.csv")




#####environmental data summaries#####
Sites <- sample_metadata %>% 
  dplyr::select(Name,Code) %>% 
  unique() %>% 
  left_join(rownames_to_column(as.data.frame.matrix(table(sample_metadata$Code, sample_metadata$Year)),"Code"))  #add in number of samples per river per year

#summary
SampleSummary <- sample_metadata %>% 
  select(Name,Code,SampleID,Year,Month,Day,Date,VolFil,MaxFlow,maxPSI,River.Width,MeanDepth,MeanFlow,WaterTemp,
         mmHg,kPa,DO.,DOmg.L,C.us..cm,pH.mV,pH,ORP.mV,Chlorine,Alkalinity,CC.) %>% 
  group_by(Year) %>% 
  summarize(MeanVolFil=mean(VolFil, na.rm = T)) ###still working on this
  
#data plots
VolHist <- ggplot(sample_metadata,aes(x=VolFil))+
  geom_histogram(binwidth=0.01,fill="grey25", color="white", alpha=0.9)+
  scale_y_continuous(limits=c(0,200),expand = expansion(mult = c(0, .02)))+
  geom_vline(xintercept=mean(sample_metadata$VolFil,na.rm=T),color="red",linetype="longdash")+
  labs(x = "Volume Filtered (L)", y = "Count", title="Liters Filtered per Sample")+
  theme_bw()+
  theme(axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(angle = 90, size=20, vjust=1),
        axis.title.x = element_text(angle = 0, size=20, vjust=-0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank());VolHist

MonthHist <- ggplot(sample_metadata,aes(x=Month,fill=Year))+
  geom_bar()+
  scale_fill_manual(values=c("dodgerblue4","steelblue3","skyblue"))+
  facet_wrap(~Year)+
  scale_y_continuous(limits=c(0,300),expand = expansion(mult = c(0, .02)))+
  labs(x = "Month", y = "Count",title="Samples Per Month")+
  theme_bw()+
  theme(axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(angle = 90, size=20, vjust=1),
        axis.title.x = element_text(angle = 0, size=20, vjust=-0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank());MonthHist







#####save workspace#####
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabarcoding_qPCR/eDNA_v_qPCR_DataProcessing_workspace.RData")
