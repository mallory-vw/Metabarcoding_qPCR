#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabarcoding_qPCR/eDNA_v_qPCR_GLMs_final_workspace.RData")

#load libraries
library("tidyverse")
library("ggeffects")
library("fitdistrplus")
library("stargazer")
library("glmmTMB")
library("DHARMa")
library("lme4")
library("bbmle") #for AICtab
library("gridExtra")
library("car")
library("MuMIn")
library("pscl")
library("brms")
library("ggmap")
library("metR")
library("tmap")
library("tmaptools")
library("RColorBrewer")
library("janitor")

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load data (see eDNA_v_qPCR_DataProcessing.R)####
AtlSalmon_data_raw <- read.csv("qPCR_eDNA_AtlSalmon_CleanedData_7Nov2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

PinkSalmon_data_raw <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_7Nov2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

ArcticCharr_data_raw <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_7Nov2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X)) 


SalmonDNAConc <- AtlSalmon_data_raw %>% 
  dplyr::select(SampleID,DNAConc) %>% 
  rename(AtlDNAConc=DNAConc) %>% 
  distinct()

PinkDNAConc <- PinkSalmon_data_raw %>% 
  dplyr::select(SampleID,DNAConc) %>% 
  rename(PinkDNAConc=DNAConc)%>% 
  distinct()

CharrDNAConc <- ArcticCharr_data_raw %>% 
  dplyr::select(SampleID,DNAConc) %>% 
  rename(CharrDNAConc=DNAConc)%>% 
  distinct()

DNAConc <- SalmonDNAConc %>% 
  full_join(PinkDNAConc, by = "SampleID") %>% 
  full_join(CharrDNAConc, by="SampleID") %>% 
  rowwise %>% 
  mutate(MatchAtlPink=identical(AtlDNAConc,PinkDNAConc),
         MatchPinkCharr=identical(PinkDNAConc, CharrDNAConc),
         MatchAtlCharr=identical(AtlDNAConc, CharrDNAConc)) 
rm(SalmonDNAConc,PinkDNAConc,CharrDNAConc)

##############All Species################
####Data Prep####
#Bind together all species data
AtlSalmon_data_int <- AtlSalmon_data_raw %>% 
  mutate(Species="Atlantic salmon")
PinkSalmon_data_int <- PinkSalmon_data_raw %>% 
  mutate(Species = "Pink salmon")
ArcticCharr_data_int <- ArcticCharr_data_raw %>% 
  mutate(Species = "Arctic charr") %>% 
  dplyr::select(!(Notes))

All_data_raw <- AtlSalmon_data_int %>% 
  rbind(PinkSalmon_data_int,ArcticCharr_data_int) %>% 
  mutate_at(c("MeanFlow","MeanDepth"), ~na_if(., 0)) %>% #replace 0s in MeanDepth and MeanFlow with NA
  mutate(RiverOutput = River.Width*MeanDepth*MeanFlow) %>% #calculate RiverOutput = Width*Depth*Flow
  relocate(RiverOutput, .after=MeanFlow) %>% 
  mutate(pHcalc = 7-(pH.mV/57.14)) %>% 
  relocate(pHcalc, .after=pH.mV)

rm(AtlSalmon_data_int,PinkSalmon_data_int,ArcticCharr_data_int)

#data summary by year
AllRiverData <- All_data_raw %>%  #first look at the River samples themselves - how many Rivers per year do we have data for?
  dplyr::select(Name,Code,Year,Date,NewLat,NewLon,#pull out only the river parameters
                "River.Width",MeanDepth,MeanFlow,RiverOutput,WaterTemp,
                pH,pH.mV,pHcalc,Alkalinity,Chlorine,
                mmHg, kPa,"DO.","DOmg.L","C.us..cm","ORP.mV","CC.") %>% 
  distinct(Name,Code,Year,Date, .keep_all = T) %>%  #remove duplicates based on River and Date
  group_by(Year)

AllRiverData %>% #generate river param summaries
  summarize(Rivers = length(unique(Name)),
            Width = sum(!is.na(River.Width)),
            Depth = sum(!is.na(MeanDepth)),
            Flow = sum(!is.na(MeanFlow)),
            Output = sum(!is.na(RiverOutput)),
            Temp = sum(!is.na(WaterTemp)),
            pH = sum(!is.na(pH)),
            Alk = sum(!is.na(Alkalinity)),
            Chlorine = sum(!is.na(Chlorine)),
            Hg = sum(!is.na(mmHg)),
            kPa = sum(!is.na(kPa)),
            "domg/L" = sum(!is.na(DOmg.L)),
            "C.us.cm" = sum(!is.na(C.us..cm)),
            ORP = sum(!is.na(ORP.mV)),
            CC = sum(!is.na(CC.))) 

#Normalize data by Vol filtered
All_data_Processing <- All_data_raw %>% 
  mutate(RawReadsPerLitre = RawReads/VolFil,
         PropRawReadsPerLitre = RawPropReads/VolFil,
         CorrectedReadsPerLitre = CorrectedReads/VolFil,
         PropCorrectedReadsPerLitre = CorrectedPropReads/VolFil,
         QuantMeanPerLitre = QuantMean/VolFil) %>% 
  relocate(RawReadsPerLitre, .after=RawReads) %>% 
  relocate(PropRawReadsPerLitre, .after = RawPropReads) %>% 
  relocate(CorrectedReadsPerLitre, .after=CorrectedReads) %>% 
  relocate(PropCorrectedReadsPerLitre, .after = CorrectedPropReads) %>% 
  relocate(QuantMeanPerLitre, .after = QuantMean) %>% 
  mutate(Code=factor(Code),
         Type=factor(Type),
         Run=factor(Run),
         Year=factor(Year),
         DNAConcScale=scale(DNAConc, center=T, scale=T)[,1]) %>% #scale the DNAConc to 0 mean and var 1
  relocate(DNAConcScale, .after = DNAConc) %>% 
  mutate(QuantMeanPerLitre=ifelse(is.na(QuantMeanPerLitre),0,QuantMeanPerLitre)) %>% #replace NA with 0 in QuantMeanPerLitre
  mutate(across(c(PropRawReadsPerLitre,PropCorrectedReadsPerLitre), ~ replace(., is.nan(.), 0))) #replace NaN with 0 in PropRawReadsPerLitre, PropCorrectedReadsPerLitre






#Samples per River before and after filtering
SamplesPerRiver <- data.frame(unclass(table(All_data_raw$Code))) %>% 
  mutate(Code=rownames(.),
         SamplesAfterFilter = data.frame(unclass(table(All_data_Processing$Code)))[,1]) %>% #Add in the Number of samples after the filter
  rename(TotalSamples = "unclass.table.All_data_raw.Code..") %>% 
  relocate(Code) #Move Code to the front


#subset out the important data for models
All_data_ModelVars <- All_data_Processing %>% 
  dplyr::select("SampleID","SampleID2","Type","Name","Code","Year","Month","Day","Date","NewLat","NewLon",
                "Marker","Taxon",
                "RawVertReadsPerSample","RawReads","RawPropReads","RawReadsPerLitre","RawPropReads","PropRawReadsPerLitre",
                "CorrectedVertReadsPerSample","CorrectedReads","CorrectedReadsPerLitre","CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFil","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result","Species",
                "River.Width","MeanDepth","MeanFlow","RiverOutput","WaterTemp","pH","pH.mV","Alkalinity","Chlorine","DOmg.L","ORP.mV")

#summary data
All_data_ModelVars %>% 
  group_by(Year) %>% 
  summarize(Samples = length(unique(SampleID))) 

#data summary by year
SampleRiverData <- All_data_ModelVars %>%  #first look at the River samples themselves - how many Rivers per year do we have data for?
  dplyr::select(Name,Code,Year,Date,NewLat,NewLon,#pull out only the river parameters
                "River.Width",MeanDepth,MeanFlow,RiverOutput,WaterTemp,
                pH,pH.mV,Alkalinity,Chlorine,
                "DOmg.L","ORP.mV") %>% 
  distinct(Name,Code,Year,Date, .keep_all = T) %>%  #remove duplicates based on River and Date
  group_by(Year)

SampleRiverData %>% #generate river param summaries
  summarize(Rivers = length(unique(Name)),
            Width = sum(!is.na(River.Width)),
            Depth = sum(!is.na(MeanDepth)),
            Flow = sum(!is.na(MeanFlow)),
            Output = sum(!is.na(RiverOutput)),
            Temp = sum(!is.na(WaterTemp)),
            pH = sum(!is.na(pH)),
            Alk = sum(!is.na(Alkalinity)),
            Chlorine = sum(!is.na(Chlorine)),
            "domg/L" = sum(!is.na(DOmg.L)),
            ORP = sum(!is.na(ORP.mV))) 

#####Site Map#####
#any sample overlap between years?
Sites_2019 <- unique(subset(All_data_ModelVars, All_data_ModelVars$Year == "2019")$Code)
Sites_2020 <- unique(subset(All_data_ModelVars, All_data_ModelVars$Year == "2020")$Code)
Sites_2021 <- unique(subset(All_data_ModelVars, All_data_ModelVars$Year == "2021")$Code)

intersect(Sites_2019,Sites_2020) #Not shared
intersect(Sites_2020,Sites_2021) #Not shared
intersect(Sites_2019,Sites_2021) #"109" "CHR" "ENG" "HUR" "IKA" "KIN" "MBB" "PAN" "SUS"

metadata_plot <- All_data_ModelVars %>% 
  distinct(Name,Code,Year,.keep_all = T) %>% 
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
       file = "eDNARivers_GeneticSiteMap_15Nov2023.pdf", 
       height = 15, 
       width = 10, 
       units = "in")



#####Examine Model data - No Outlier Removal#####
All_data_ModelVars_Prop <- All_data_ModelVars %>%  #pull out just the Proportion data for modelling
  dplyr::select(!c(RawVertReadsPerSample,RawReads,RawPropReads,RawReadsPerLitre,
                   CorrectedVertReadsPerSample,CorrectedReads,CorrectedReadsPerLitre,CorrectedPropReads)) 

#Resonse Variable: Prop Raw reads per liter
hist(All_data_ModelVars_Prop$PropRawReadsPerLitre) #data is not normal
format(range(All_data_ModelVars_Prop$PropRawReadsPerLitre),scientific=F) #[1] "0.000000" "1.200385", proportion can't be higher than 1
boxplot(All_data_ModelVars_Prop$PropRawReadsPerLitre, horizontal = T)

#have a proportion over 1, so remove that value
All_data_ModelVars_RawProp <- All_data_ModelVars_Prop %>% 
  filter(PropRawReadsPerLitre < 1) #removed 11 rows

#explanatory variable - Mean Quant score
hist(All_data_ModelVars_RawProp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(All_data_ModelVars_RawProp$QuantMeanPerLitre, na.rm=T) ##[1]   0.0000 345.7539
boxplot(All_data_ModelVars_RawProp$QuantMeanPerLitre, horizontal = T)

#export data
write.csv(All_data_ModelVars_RawProp,
          "AllSpecies_AllMarkers_RawProp_15Nov2023.csv")


#Resonse Variable: Prop Corrected reads per liter
hist(All_data_ModelVars_Prop$PropCorrectedReadsPerLitre) #data is not normal
format(range(All_data_ModelVars_Prop$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.000000" "1.315789", proportion can't be higher than 1
boxplot(All_data_ModelVars_Prop$PropCorrectedReadsPerLitre, horizontal = T)

#have a proportion over 1, so remove that value
All_data_ModelVars_CorrProp <- All_data_ModelVars_Prop %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 19 rows

#explanatory variable - Mean Quant score
hist(All_data_ModelVars_CorrProp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(All_data_ModelVars_CorrProp$QuantMeanPerLitre, na.rm=T) ##[1]   0.0000 345.7539
boxplot(All_data_ModelVars_CorrProp$QuantMeanPerLitre, horizontal = T)

#export data
write.csv(All_data_ModelVars_CorrProp,
          "AllSpecies_AllMarkers_CorrectedProp_15Nov2023.csv")


######Extract all necessary Vars and scale to 0 mean and Unit variance #####
#raw proportions
All_data_ModelVars_RawProp_Model <- All_data_ModelVars_RawProp %>% 
  mutate(QuantMeanPerLiterScale=scale(QuantMeanPerLitre, center=T, scale=T)[,1],
         RiverOutputScale=scale(RiverOutput, center=T, scale=T)[,1],
         WaterTempScale=scale(WaterTemp, center=T, scale=T)[,1],
         pHScale=scale(pH, center=T, scale=T)[,1],
         AlkalinityScale=scale(Alkalinity, center=T, scale=T)[,1],
         ChlorineScale=scale(Chlorine, center=T, scale=T)[,1],
         DOScale=scale(DOmg.L, center=T,scale=T)[,1],
         ORPScale=scale(ORP.mV,center=T,scale=T)[,1]) %>% 
  dplyr::select(SampleID, PropRawReadsPerLitre,QuantMeanPerLiterScale,DNAConcScale,WaterTempScale,DOScale,ORPScale,
                Taxon,Year,Type,Marker,Result,Run,Code) %>% 
  filter(complete.cases(.)) #2735 data points

#corrected proportions
All_data_ModelVars_CorrProp_Model <- All_data_ModelVars_CorrProp %>% 
  mutate(QuantMeanPerLiterScale=scale(QuantMeanPerLitre, center=T, scale=T)[,1],
         RiverOutputScale=scale(RiverOutput, center=T, scale=T)[,1],
         WaterTempScale=scale(WaterTemp, center=T, scale=T)[,1],
         pHScale=scale(pH, center=T, scale=T)[,1],
         AlkalinityScale=scale(Alkalinity, center=T, scale=T)[,1],
         ChlorineScale=scale(Chlorine, center=T, scale=T)[,1],
         DOScale=scale(DOmg.L, center=T,scale=T)[,1],
         ORPScale=scale(ORP.mV,center=T,scale=T)[,1]) %>% 
  dplyr::select(SampleID, PropCorrectedReadsPerLitre,QuantMeanPerLiterScale,DNAConcScale,WaterTempScale,DOScale,ORPScale,
                Taxon,Year,Type,Marker,Result,Run,Code) %>% 
  filter(complete.cases(.)) #2727 data points


#####All Markers Variable Selection - Raw#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
Full_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                        PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Run) + (1|Code),
                        na.action = na.fail,
                        family = beta_family(),
                        ziformula = ~1, #ziformula=~1 means probability of structural zero is the same for all obs
                        REML=T)


summary(Full_LMM_Raw)
AICc(Full_LMM_Raw) #-3040.086



#sort out random using AIC
Random1_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                              PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result,
                           na.action = na.fail,
                           family = beta_family(),
                              ziformula = ~1, #ziformula=~1 means probability of structural zero is the same for all obs
                              REML=T)
Random2_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                              PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Run),
                           na.action = na.fail,
                           family = beta_family(),
                              ziformula = ~1, #ziformula=~1 means probability of structural zero is the same for all obs
                              REML=T)
Random3_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                              PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Code),
                           na.action = na.fail,
                           family = beta_family(),
                              ziformula = ~1, #ziformula=~1 means probability of structural zero is the same for all obs
                              REML=T)

#check AICc of the models
AICc(Full_LMM_Raw)  #-3040.086  
AICc(Random1_LMM_Raw) #-2934.589
AICc(Random2_LMM_Raw)  #  -2966.375
AICc(Random3_LMM_Raw) # -3022.53

#Full Model has the best AIC score

#zero-inflation structure
Zi1_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                       PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Run) + (1|Code),
                       na.action = na.fail,
                       family = beta_family(),
                       ziformula = ~Result, #probability of structural zero varies by Result
                       REML=T)
Zi2_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                       PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Run) + (1|Code),
                       na.action = na.fail,
                       family = beta_family(),
                       ziformula = ~DNAConcScale, #probability of structural zero varies by DNAConcScale
                       REML=T)
Zi3_LMM_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                       PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Run) + (1|Code),
                       na.action = na.fail,
                       family = beta_family(),
                       ziformula = ~., #probability of structural zero is the same as the main formula
                       REML=T)

#check AICc of the models
AICc(Full_LMM_Raw)  #-3040.086  
AICc(Zi1_LMM_Raw) #-4522.138
AICc(Zi2_LMM_Raw)  #  -3039.257
AICc(Zi3_LMM_Raw) # -5466.005

#Zi3_LMM_Raw has the best AIC score

#Select Fixed Effects, use AIC and DO NOT use REML
Full_LMM_Raw_Fixed <-  glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                               PropRawReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + WaterTempScale + DOScale + ORPScale + Taxon + Year + Type + Marker + Result + (1|Run) + (1|Code),
                               na.action = na.fail,
                               family = beta_family(),
                               ziformula = ~., #probability of structural zero is the same as the main formula
                               REML=F)
AICc(Full_LMM_Raw_Fixed) #-5514.874

#use dredge to run all possible sub-models
Dredge_Full_LMM_Raw_Fixed <- dredge(Full_LMM_Raw_Fixed)
print(Dredge_Full_LMM_Raw_Fixed)
# Model selection table 
#      cnd((Int)) zi((Int)) dsp((Int)) cnd(DNA)   cnd(DOS) cnd(Mrk)   cnd(ORP) cnd(QMP) cnd(Rsl) cnd(Txn) cnd(Typ)   cnd(WTS) cnd(Yer) df   logLik    AICc  delta weight
# 125     -1.2130    0.2092          +                            +  0.0663100   0.2912        +        +                              15 1565.245 -3100.3   0.00  0.098
# 117     -1.2060    0.2092          +                            +              0.2896        +        +                              14 1564.096 -3100.0   0.27  0.086
# 118     -1.2130    0.2092          + -0.07101                   +              0.2925        +        +                              15 1565.020 -3099.9   0.45  0.079
# 126     -1.2190    0.2092          + -0.06457                   +  0.0611000   0.2937        +        +                              16 1566.000 -3099.8   0.51  0.076
# 127     -1.2180    0.2092          +           0.0566000        +  0.0797900   0.2924        +        +                              16 1565.617 -3099.0   1.28  0.052
# 381     -1.2110    0.2092          +                            +  0.0665400   0.2908        +        +           1.115e-02          16 1565.265 -3098.3   1.98  0.037

#Top 6 models are basically equivalent, based on deltaAIC (<2)
#use the simplest model, with just QuantMeanPerLitre,Marker,Result,Taxon


##### AllMarkers Final model#####
FinalModel_Raw <- glmmTMB(data = All_data_ModelVars_RawProp_Model, 
                          PropRawReadsPerLitre ~ QuantMeanPerLiterScale + Taxon + Marker + Result + (1|Run) + (1|Code),
                          na.action = na.fail,
                          family = beta_family(),
                          ziformula = ~1, #ziformula=~1 means probability of structural zero is the same for all obs
                          REML=T)
plot(simulateResiduals(FinalModel_Raw))
summary(FinalModel_Raw) 
Anova(FinalModel_Raw,type="III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: PropCorrectedReadsPerLitre
#                           Chisq Df Pr(>Chisq)    
#   (Intercept)            14.398  1  0.0001480 ***
#   QuantMeanPerLiterScale 64.976  1  7.581e-16 ***
#   DNAConcScale           12.617  1  0.0003822 ***
#   Marker                 22.347  2  1.404e-05 ***
#   Result                 28.509  3  2.840e-06 ***
#   Taxon                   6.329  2  0.0422345 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#plain LM to compare
LM_finalmodel_NoZero <- lm(data = All_data_NoZero_Prop_NoOutliers_Model, 
                           na.action = na.omit,
                           PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + Marker + Result + Taxon + Run + Code)
plot(simulateResiduals(LM_finalmodel_NoZero))
summary(LM_finalmodel_NoZero) 

#Final model with log normal
FinalModel_log_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model,
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + Marker + Result + Taxon + (1|Run) + (1|Code),
                                 na.action = na.omit,
                                 REML=T,
                                 family = gaussian(link="logit"))
plot(simulateResiduals(FinalModel_log_NoZero))
summary(FinalModel_log_NoZero) 


AICctab(FinalModel_NoZero,
        LM_finalmodel_NoZero,
        FinalModel_log_NoZero,
        delta=T,base=T)
#                       AICc   dAICc  df 
# FinalModel_NoZero     -816.4    0.0 13 
# LM_finalmodel_NoZero  -313.9  502.5 141
# FinalModel_log_NoZero -287.1  529.3 13 

#######Split by species#######


















#####plot the models####
# Extract the prediction data frame
pred_mm_NoZero <- ggpredict(FinalModel_NoZero, terms = c("QuantMeanPerLiterScale"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_NoZero_plot <- ggplot(pred_mm_NoZero) + 
  geom_line(data = pred_mm_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_mm_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_NoZero_Prop_NoOutliers_Model, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLiterScale, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "All Species Beta GLMM\nPropCorrectedReadsPerLitre ~ scale(QuantMeanPerLitre) + scale(DNAConc) + Marker + Species + Result + (1|qPCRRun) + (1|River)") +
  theme_bw() +
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalModel_NoZero_plot

#export plots
ggsave(FinalModel_NoZero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_NoZero_NoOutliers_GLMMBeta_14Nov2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width









######save workspace######
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabarcoding_qPCR/eDNA_v_qPCR_GLMs_final_workspace.RData")
