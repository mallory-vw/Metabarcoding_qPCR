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

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load data (see eDNA_v_qPCR_DataProcessing.R)####
AtlSalmon_data_raw <- read.csv("qPCR_eDNA_AtlSalmon_CleanedData_7Nov2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

PinkSalmon_data_raw <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_7Nov2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

ArcticCharr_data_raw <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_7Nov2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X)) 

###############Models with 0 read count removed###########
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
  rbind(PinkSalmon_data_int,ArcticCharr_data_int)
rm(AtlSalmon_data_int,PinkSalmon_data_int,ArcticCharr_data_int)

#calculate pH from pH.mV
#pH = 7-(mV/X)
#need to calculate the constant in the formula for each river
All_data_pHcalc <- All_data_raw %>% 
  dplyr::select("Name","Code","Year","Month","Day","Date",
                "NewLat","NewLon","WaterTemp","pH.mV","pH") %>% 
  arrange(Date) %>% 
  group_by(Name,Date) %>% 
  distinct(Name,Date, .keep_all = T) %>% 
  mutate(pHConst=pH.mV/(7-pH)) %>% #calculate the constant for each sample that has both pH and pH.mV
  mutate(pHConst2=case_when(Year == 2020 ~ mean(filter(., Year==2019)$pHConst)),  #set the constant for 2020 to be the mean constant of 2019
         pH2=case_when(Year == 2020 ~ 7-(pH.mV/pHConst2)), #calculate pH for 2020 samples using the average constant from 2019
         pHfinal=coalesce(pH,pH2)) %>% 
  relocate(pHfinal, .after= WaterTemp)

#look at ranges for the measured and calculated pH values
All_data_pHcalc %>% 
  group_by(Year) %>% 
  summarize(pHmin = min(pHfinal),
            pHmax = max(pHfinal))

All_data_pHcalc %>% 
  filter(pHfinal >= 7.7)

#I'm not sure I trust the calculations for 2020 - the highest pH is 8.42, which seems SO HIGH






#Normalize data by Vol filtered and remove any rows with 0 corrected reads
All_data_NoZero <- All_data_raw %>% 
  filter(CorrectedReads > 0) %>% 
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
         DNAConcScale=scale(DNAConc, center=T, scale=T)[,1]) %>% 
  relocate(DNAConcScale, .after = DNAConc)

#subset out the important data for models
All_data_NoZero_Model <- All_data_NoZero %>% 
  dplyr::select("SampleID","SampleID2","Type","Name","Code","Year","Month","Day","Date","NewLat","NewLon",
                "Marker","Taxon",
                "RawVertReadsPerSample","RawReads","RawPropReads","RawReadsPerLitre","RawPropReads","PropRawReadsPerLitre",
                "CorrectedVertReadsPerSample","CorrectedReads","CorrectedReadsPerLitre","CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFil","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result","Species",
                "River.Width","MeanDepth","MeanFlow","WaterTemp","pH","pH.mV","Alkalinity","Chlorine") %>% 
  mutate_at(c("MeanFlow","MeanDepth"), ~na_if(., 0)) %>% #replace 0s in MeanDepth and MeanFlow with NA
  mutate(RiverOutput = River.Width*MeanDepth*MeanFlow) %>% #calculate RiverOutput = Width*Depth*Flow
  relocate(RiverOutput, .after=MeanFlow) %>% 
  mutate(pHcalc = 7-(pH.mV/57.14)) %>% 
  relocate(pHcalc, .after=pH.mV)



#####remove outliers#####
All_data_NoZero_Prop <- All_data_NoZero_Model %>%  #pull out just the Proportion data for modelling
  dplyr::select(!c(RawVertReadsPerSample,RawReads,RawPropReads,RawReadsPerLitre,
                   CorrectedVertReadsPerSample,CorrectedReads,CorrectedReadsPerLitre,CorrectedPropReads)) 

#Resonse Variable: Prop corrected reads per liter
hist(All_data_NoZero_Prop$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score#Response variable - proportion of total reads
format(range(All_data_NoZero_Prop$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.0000007896875" "1.3157894736842"
boxplot(All_data_NoZero_Prop$PropCorrectedReadsPerLitre)
summary(All_data_NoZero_Prop$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.0815611 0.2431401 0.3361202 0.5313641 1.3157895 
IQR <- IQR(All_data_NoZero_Prop$PropCorrectedReadsPerLitre)
quartiles <- quantile(All_data_NoZero_Prop$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_Prop_NoOutliers_Temp <- All_data_NoZero_Prop %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper) %>%  #2 outliers removed
  dplyr::select(!PropRawReadsPerLitre)

rm(IQR, quartiles,Lower,Upper)
hist(All_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(All_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre)

#explanatory variable - Mean Quant score
hist(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T) ##Prop corrected reads per liter
boxplot(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)
summary(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.6029   6.3435  12.1762  18.8963  21.1311 248.6284      117 
IQR <- IQR(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_Prop_NoOutliers <- All_data_NoZero_Prop_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper) #182 outliers removed
rm(IQR, quartiles,Lower,Upper,All_data_NoZero_Prop_NoOutliers_Temp)
hist(All_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre) 
boxplot(All_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre)

#still have a proportion over 1, so remove that value
format(range(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F)
All_data_NoZero_Prop_NoOutliers <- All_data_NoZero_Prop_NoOutliers %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 14 rows

#export data
write.csv(All_data_NoZero_Prop_NoOutliers,
          "AllSpecies_AllMarkers_NoZero_NoOutliers_14Nov2023.csv")

#Resonse Variable: Prop raw reads per liter
hist(All_data_NoZero_Prop$PropRawReadsPerLitre) #data is not normal, can't use z-score#Response variable - proportion of total reads
format(range(All_data_NoZero_Prop$PropRawReadsPerLitre),scientific=F) #[1] [1] "0.0000007895833" "1.2003846153846"
boxplot(All_data_NoZero_Prop$PropRawReadsPerLitre)
summary(All_data_NoZero_Prop$PropRawReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.0833300 0.2633605 0.3272967 0.5016354 1.2003846 
IQR <- IQR(All_data_NoZero_Prop$PropRawReadsPerLitre)
quartiles <- quantile(All_data_NoZero_Prop$PropRawReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_RawProp_NoOutliers_Temp <- All_data_NoZero_Prop %>% 
  filter(PropRawReadsPerLitre > Lower & PropRawReadsPerLitre < Upper) %>%  #1 outlier removed
  dplyr::select(!PropCorrectedReadsPerLitre)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_NoZero_RawProp_NoOutliers_Temp$PropRawReadsPerLitre) 
boxplot(All_data_NoZero_RawProp_NoOutliers_Temp$PropRawReadsPerLitre)

#explanatory variable - Mean Quant score
hist(All_data_NoZero_RawProp_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(All_data_NoZero_RawProp_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T) #[1]   0.6029412 248.6284000
boxplot(All_data_NoZero_RawProp_NoOutliers_Temp$QuantMeanPerLitre)
summary(All_data_NoZero_RawProp_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.6029   6.3626  12.1762  18.8943  21.1311 248.6284      117 
IQR <- IQR(All_data_NoZero_RawProp_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(All_data_NoZero_RawProp_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_RawProp_NoOutliers <- All_data_NoZero_RawProp_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper) #182 outliers removed
rm(IQR, quartiles,Lower,Upper,All_data_NoZero_RawProp_NoOutliers_Temp)
hist(All_data_NoZero_RawProp_NoOutliers$QuantMeanPerLitre) 
boxplot(All_data_NoZero_RawProp_NoOutliers$QuantMeanPerLitre)

#still have a proportion over 1, so remove that value
format(range(All_data_NoZero_RawProp_NoOutliers$PropRawReadsPerLitre),scientific=F)
All_data_NoZero_RawProp_NoOutliers <- All_data_NoZero_RawProp_NoOutliers %>% 
  filter(PropRawReadsPerLitre < 1) #removed 7 rows

#export data
write.csv(All_data_NoZero_RawProp_NoOutliers,
          "AllSpecies_AllMarkers_NoZero_Raw_NoOutliers_14Nov2023.csv")

######Extract all necessary Vars and scale to 0 mean and Unit variance #####
All_data_NoZero_Prop_NoOutliers_Model <- All_data_NoZero_Prop_NoOutliers %>% 
  mutate(#PropCorrectedReadsPerLitreScale=scale(PropCorrectedReadsPerLitre, center=T, scale=T)[,1],
         QuantMeanPerLiterScale=scale(QuantMeanPerLitre, center=T, scale=T)[,1],
         RiverOutputScale=scale(RiverOutput, center=T, scale=T)[,1],
         WaterTempScale=scale(WaterTemp, center=T, scale=T)[,1],
         pHScale=scale(pH, center=T, scale=T)[,1],
         AlkalinityScale=scale(Alkalinity, center=T, scale=T)[,1],
         ChlorineScale=scale(Chlorine, center=T, scale=T)[,1]) %>% 
  dplyr::select("SampleID","Type","Code","Year","Month","Day","Date","NewLat","NewLon","Marker","Taxon",
                "PropCorrectedReadsPerLitre","DNAConcScale","QuantMeanPerLiterScale","Result","Run",
                "RiverOutputScale","WaterTempScale","pHScale","AlkalinityScale","ChlorineScale")

#####All Markers Variable Selection#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
Full_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model, 
                           PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + Marker + Type + Result + Taxon + DNAConcScale + WaterTempScale + Year + (1|Run) + (1|Code),
                           na.action = na.omit,
                           family = beta_family(),
                           REML=T)
summary(Full_LMM_NoZero)
AICc(Full_LMM_NoZero) #-797.9297

#sort out random using AIC
Random1_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model, 
                              PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + Marker + Type + Result + Taxon + DNAConcScale + WaterTempScale + Year,
                              na.action = na.omit,
                              family = beta_family(),
                              REML=T)
Random2_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model, 
                              PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + Marker + Type + Result + Taxon + DNAConcScale + WaterTempScale + Year + (1|Run),
                              na.action = na.omit,
                              family = beta_family(),
                              REML=T)
Random3_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model, 
                              PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + Marker + Type + Result + Taxon + DNAConcScale + WaterTempScale + Year + (1|Code),
                              na.action = na.omit,
                              family = beta_family(),
                              REML=T)

#check AICc of the models
AICc(Full_LMM_NoZero)  #-797.9297  
AICc(Random1_LMM_NoZero) #-651.4994
AICc(Random2_LMM_NoZero)  #  -669.2169
AICc(Random3_LMM_NoZero) # -788.4696

#Full Model has the best AIC score

#Select Fixed Effects, use AIC and DO NOT use REML
Full_LMM_NoZero_Fixed <-  glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model, 
                                  PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + Marker + Type + Result + Taxon + DNAConcScale + WaterTempScale + Year + (1|Run) + (1|Code),
                                  na.action = na.fail,
                                  family = beta_family(),
                                  REML=F)
AICc(Full_LMM_NoZero_Fixed) #-828.4839

#use dredge to run all possible sub-models
Dredge_Full_LMM_NoZero_Fixed <- dredge(Full_LMM_NoZero_Fixed)
print(Dredge_Full_LMM_NoZero_Fixed)
# Model selection table 
#     cnd((Int)) dsp((Int)) cnd(DNA) cnd(Mrk) cnd(QMP) cnd(Rsl) cnd(Txn) cnd(Typ)   cnd(WTS) cnd(Yer) df  logLik   AICc  delta weight
# 32     -0.5809          + -0.24580        +   0.4399        +        +                              13 431.831 -837.2   0.00  0.267
# 160    -0.8544          + -0.26290        +   0.4486        +        +                            + 15 433.690 -836.8   0.43  0.215
# 96     -0.5671          + -0.25080        +   0.4408        +        +           0.0761400          14 432.257 -836.0   1.22  0.145
# 144    -0.7418          + -0.27060        +   0.4429        +                                     + 13 430.751 -835.0   2.16  0.091

#Top 3 models are basically equivalent, based on deltaAIC (<2)
#use the simplest model, with just QuantMeanPerLitre,DNA Concentration,Marker,Result,Taxon


##### AllMarkers Final model#####
FinalModel_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers_Model, 
                             PropCorrectedReadsPerLitre ~ QuantMeanPerLiterScale + DNAConcScale + Marker + Result + Taxon + (1|Run) + (1|Code),
                             na.action = na.omit,
                             family = beta_family(),
                             REML=T)
plot(simulateResiduals(FinalModel_NoZero))
summary(FinalModel_NoZero) 
Anova(FinalModel_NoZero,type="III") 
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
