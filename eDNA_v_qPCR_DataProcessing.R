#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabarcoding_qPCR/eDNA_v_qPCR_DataProcessing_workspace.RData")

#load libraries
library("tidyverse")
library("ggmap")
library("metR")
library("tmap")
library("tmaptools")
library("RColorBrewer")

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load all raw data#####
metadata <- read.csv("RiverSamplingMetadata_2019_2021.csv",na.strings=c("","NA")) %>% 
  rename(Name=RiverName,
         Code=RiverCode) %>% 
  mutate(Year=factor(Year))

sample_metadata_raw <- read.csv("SampleMetadata_EnviroData_7Nov2023.csv",na.strings=c("","NA")) %>% #NOTE: Date edited in Excel to Month/Day/Year
  rename(Name=RiverName,
         Code=RiverCode,
         Verified=Sample.Verified) %>% 
  mutate(Year=factor(Year)) 
  
site_metadata_raw <- read.csv("SiteMetadata_EnviroData_7Nov2023.csv",na.strings=c("","NA")) %>% #NOTE: Date as Month/Day/Year
  rename(Name=RiverName,
         Code=RiverCode) %>% 
  mutate(Year=factor(Year)) 


eDNA_data_raw <- read.csv("eDNA_results_PrelimCleaned_5Oct2023.csv", na.strings = "N/A") %>% 
  mutate(CorrectedDepth = as.integer(CorrectedDepth),
         Type=factor(Type)) %>% 
  rename(RawVertReadsPerSample=RawVertDepthPerSample,
         RawReads=RawDepth,
         CorrectedReads=CorrectedDepth,
         SampleID=Sample)

AtlSalmon_qPCR_raw <- read.csv("qPCR_SalmoSalar_results_PrelimCleaned_5Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID,
         DNAConc=DNAConc_pg_uL)
PinkSalmon_qPCR_raw <- read.csv("qPCR_PinkSalmon_results_PrelimCleaned_5Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID,
         DNAConc=DNAConc_pg_uL)
ArcticCharr_qPCR_raw <- read.csv("qPCR_ArcticCharr_results_PrelimCleaned_10Oct2023.csv") %>% 
  rename(SampleID=CEGA_Sample_ID,
         DNAConc=DNAConc_pg_uL)

#####Organize Sample Site data#####
sample_metadata <- sample_metadata_raw %>% 
  left_join(site_metadata_raw, by=c("Code","Year","Date")) %>%  #merge in site data, manually changed Crooked River, Shinney's River side channel codes to match properly
  rename(SampleName = Name.x, #rename fields
         SiteName = Name.y,
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
         Month=factor(Month, levels=c("August","September","October","November","December"))) %>% 
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
  dplyr::select(!c(Year2,NewLatPrep,NewLonPrep,SampleLat,SampleLon,SiteLat,SiteLon,
                   Left.Latitude,Left.Longitude,Center.Latitude,Center.Longitude,Right.Latitude,Right.Longitude)) %>%  #remove extraneous columns
  relocate(c(Month,Day),.after=Year) %>% 
  relocate(SiteName, .after = SampleName) %>%
  relocate(c(NewLat,NewLon), .before=WPT) %>% 
  relocate(c(MeanDepth,MeanFlow), .after=River.Width) %>% 
  relocate(c(SampleNotes,SiteNotes), .before=Depth1) %>% 
  relocate(FullTrans, .after=WC)
  



#####Site Map#####
#any sample overlap between years?
Sites_2019 <- unique(subset(sample_metadata, sample_metadata$Year == "2019")$Code)
Sites_2020 <- unique(subset(sample_metadata, sample_metadata$Year == "2020")$Code)
Sites_2021 <- unique(subset(sample_metadata, sample_metadata$Year == "2021")$Code)

intersect(Sites_2019,Sites_2020) #Not shared
intersect(Sites_2020,Sites_2021) #Not shared
intersect(Sites_2019,Sites_2021) #"ENG" "FOR" "HUR" "IKA" "KIN" "MBB" "PAN" "PIN" "109" "SAN" "CHR" "SUS"

metadata_plot <- sample_metadata %>% 
  distinct(SampleName,Code,Year,.keep_all = T) %>% 
  group_by(Code) %>%
  mutate(CodeYear=ifelse(anyDuplicated(Code),"2019 and 2021",as.character(Year))) %>% #if any sites have both 2019 and 2020, note them
  distinct(Code,CodeYear,.keep_all = T) %>% 
  mutate(CodeYear=factor(CodeYear, levels = c("2019","2020","2021","2019 and 2021"))) %>% 
  drop_na(NewLat)


#Create a bounding box using tmaps
box <- bb(x=c(min(metadata_plot$NewLon,na.rm=T), #Lonmin -64.1 
             min(metadata_plot$NewLat,na.rm=T), #Latmin 46.7
             max(metadata_plot$NewLon,na.rm=T), #Lonmax -52.8
             max(metadata_plot$NewLat,na.rm=T))) #Latmax 59.4
lat_centre <- ((max(metadata_plot$NewLat,na.rm=T) - min(metadata_plot$NewLat,na.rm=T))/2)+min(metadata_plot$NewLat,na.rm=T)
lon_centre <- ((max(metadata_plot$NewLon,na.rm=T) - min(metadata_plot$NewLon,na.rm=T))/2)+min(metadata_plot$NewLon,na.rm=T)

#get googlemap
base <- get_googlemap(center = c(lon_centre,lat_centre), 
                      zoom = 5,scale = 2, #3 = continent, 21=building
                      maptype = "satellite",
                      language = "en-EN", sensor = FALSE, messaging = FALSE, 
                      urlonly = FALSE, filename = "ggmapTemp",
                      color = "color",
                      style="feature:administrative|visibility:off&style=feature:road|visibility:off")

eDNARivers_SiteMap <- ggmap(base) +
  geom_point(data=metadata_plot,size=4,#position=position_jitter(width = 0.1, height = 0.2),
             aes(x=NewLon,y=NewLat,fill=CodeYear,shape=CodeYear))+
  scale_shape_manual(values=c(21,22,23,24))+
  scale_fill_manual(values=c("#BC3C29","#4DBBD5","#E18727","#42B540"))+

  scale_x_continuous(breaks=c(-65,-60,-55), limits=c(-65,-52), labels =c(expression(paste("-65",degree)),
                                                                         expression(paste("-60",degree)),
                                                                         expression(paste("-55",degree)))) +
  scale_y_continuous(breaks=c(50,55,60),limits=c(46,59.9),labels =c(expression(paste("50",degree)),
                                                                    expression(paste("55",degree)),
                                                                    expression(paste("60",degree))))+  
  labs(x = "Longitude", y = "Latitude",fill="Sampling Year",shape="Sampling Year")+
  annotate(geom = "text", x = -62.5, y = 52.5, size = 8, label = "Labrador", colour="grey70")+
  annotate(geom = "text", x = -56, y = 48.7, size = 8, label = "NL", colour="grey70")+
  theme_bw()+
  theme(legend.position = c(0.75,0.9), 
        legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",colour="black"),
        legend.title = element_text(size=24), 
        legend.text = element_text(size=20),
        legend.key = element_blank(),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(angle = 90, size=20, vjust=1),
        axis.title.x = element_text(angle = 0, size=20, vjust=-0.5),
        panel.border = element_rect(linewidth =0.5),
        plot.margin = unit(c(5,5,7,5), "mm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(eDNARivers_SiteMap, 
       file = "eDNARivers_SiteMap_7Nov2023.pdf", 
       height = 15, 
       width = 10, 
       units = "in")

######eDNA#####
#organize total reads per sample (reads of ALL vertebrates in the sample)
eDNA_ReadsSample <- eDNA_data_raw %>% 
  select(SampleID, Type, Marker, RawVertReadsPerSample) %>% 
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
  mutate(CorrectedPropReads = as.character(signif((CorrectedReads/RawVertReadsPerSample), digits=4))) %>% 
  mutate(RawPropReads = as.character(signif((RawReads/RawVertReadsPerSample), digits=4)))

#separate data by species
AtlSalmon_eDNA_data <- eDNA_data %>% 
  filter(Taxon == "Salmo salar")

PinkSalmon_eDNA_data <- eDNA_data %>% 
  filter(Taxon == "Oncorhynchus gorbuscha")

Charr_eDNA_data <- eDNA_data %>% 
  filter(Taxon == "Salvelinus alpinus")


#####qPCR and eDNA#####
#Atlantic salmon#
#check and remove any rows that had inhibition
unique(AtlSalmon_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

AtlSalmon_qPCR_data <- AtlSalmon_qPCR_raw %>% 
  select(!IPCResult) %>% #remove the IPCResults column
  mutate(QuantMean=rowMeans(select(., c("Quant1","Quant2","Quant3")), na.rm = TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)


#join data, keeping ONLY the samples in both qPCR and eDNA results
AtlSalmon_data <- AtlSalmon_eDNA_data %>% 
  inner_join(AtlSalmon_qPCR_data) %>% 
  left_join(sample_metadata) %>% 
  relocate(VolFil, .before=Sample_Year) %>% 
  relocate(RawPropReads, .after=RawReads) %>% 
  relocate(CorrectedPropReads, .after = CorrectedReads) %>% 
  relocate(c(SampleName,Code,ReplicateSet,Year,Month,Day,Date,NewLat,NewLon),.after=Type) %>% 
  select(!Sample_Year)


#export
write.csv(AtlSalmon_data, "qPCR_eDNA_AtlSalmon_CleanedData_7Nov2023.csv")

#Pink salmon#
#check and remove any rows that had inhibition
unique(PinkSalmon_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

PinkSalmon_qPCR_data <- PinkSalmon_qPCR_raw %>% 
  select(!IPCResult) %>% #remove the IPCResults column
  mutate(QuantMean=rowMeans(select(., c("Quant1","Quant2","Quant3")), na.rm = TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)

#join data, keeping ONLY the samples in both qPCR and eDNA results
PinkSalmon_data <- PinkSalmon_eDNA_data %>% 
  inner_join(PinkSalmon_qPCR_data) %>% 
  left_join(sample_metadata) %>% 
  relocate(VolFil, .before=Sample_Year) %>% 
  relocate(RawPropReads, .after=RawReads) %>% 
  relocate(CorrectedPropReads, .after = CorrectedReads) %>% 
  relocate(c(SampleName,Code,ReplicateSet,Year,Month,Day,Date,NewLat,NewLon),.after=Type) %>% 
  select(!Sample_Year)

#export
write.csv(PinkSalmon_data, "qPCR_eDNA_PinkSalmon_CleanedData_7Nov2023.csv")

#Arctic charr#
#check and remove any rows that had inhibition
unique(ArcticCharr_qPCR_raw$IPCResult) #Make sure all samples had no inhibition

ArcticCharr_qPCR_data <- ArcticCharr_qPCR_raw %>% 
  select(!IPCResult) %>% #remove the IPCResults column
  mutate(QuantMean=rowMeans(select(., c("Quant1","Quant2","Quant3")), na.rm = TRUE)) %>% 
  relocate(QuantMean, .after=Quant3) %>% 
  filter(!SampleID %in% c("GLPR","SWA2019C2","MBB2019C2")) %>%  #remove problem samples based on the Notes field
  mutate(Ct2 = replace(Ct2, which(SampleID %in% c("TOM2019C3","PAN2019C1","PAN2019C3","SWA2019C3")), NA)) %>% #replace empty replicates with NA based on Notes field
  mutate(Ct3 = replace(Ct3, which(SampleID %in% c("TOM2019C3","PAN2019C1","PAN2019C3","SWA2019C3","SWA2019C1","KIN2019C1")), NA)) #replace empty replicates with NA based on Notes field

#join data, keeping ONLY the samples in both qPCR and eDNA results
ArcticCharr_data <- Charr_eDNA_data %>% 
  inner_join(ArcticCharr_qPCR_data) %>% 
  left_join(sample_metadata) %>% 
  relocate(VolFil, .before=Sample_Year) %>% 
  relocate(RawPropReads, .after=RawReads) %>% 
  relocate(CorrectedPropReads, .after = CorrectedReads) %>% 
  relocate(c(SampleName,Code,ReplicateSet,Year,Month,Day,Date,NewLat,NewLon),.after=Type) %>% 
  select(!Sample_Year)



#export
write.csv(ArcticCharr_data, "qPCR_eDNA_ArcticCharr_CleanedData_7Nov2023.csv")


#####save workspace#####
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabarcoding_qPCR/eDNA_v_qPCR_DataProcessing_workspace.RData")
