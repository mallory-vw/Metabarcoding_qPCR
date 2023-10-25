#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")

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
AtlSalmon_data_raw <- read.csv("qPCR_eDNA_SalmoSalar_CleanedData_19Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X)) %>% 
  rename(TotalVertReadsPerSample=RawReadsPerSample)

PinkSalmon_data_raw <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_19Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X)) %>% 
  rename(TotalVertReadsPerSample=RawReadsPerSample)

ArcticCharr_data_raw <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_19Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X)) %>% 
  rename(TotalVertReadsPerSample=RawReadsPerSample)


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

#Normalize data by Vol filtered and remove any rows with 0 corrected reads
All_data_NoZero <- All_data_raw %>% 
  filter(CorrectedReads > 0) %>% 
  mutate(RawReadsPerLitre = RawReads/VolFiltered,
         PropRawReadsPerLitre = RawPropReads/VolFiltered,
         CorrectedReadsPerLitre = CorrectedReads/VolFiltered,
         PropCorrectedReadsPerLitre = CorrectedPropReads/VolFiltered,
         QuantMeanPerLitre = QuantMean/VolFiltered) %>% 
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
  dplyr::select("SampleID","Type","Name","Code","Date","Marker","Taxon",
                "TotalVertReadsPerSample","RawReads","RawReadsPerLitre",
                "RawPropReads","PropRawReadsPerLitre",
                "CorrectedReads","CorrectedReadsPerLitre",
                "CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFiltered","DNAConc","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result","Species")



#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6


P1_NoZero <- ggplot(All_data_NoZero_Model,
             aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P2_NoZero <- ggplot(All_data_NoZero_Model,
             aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Corrected Reads Per Litre Filtered",
       title="P2")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero <- ggplot(All_data_NoZero_Model,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())#changes the axes, etc

P4_NoZero <- ggplot(All_data_NoZero_Model,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Corrected Reads Per Litre Filtered",
       title="P4") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero <- ggplot(All_data_NoZero_Model,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P6_NoZero <- ggplot(All_data_NoZero_Model,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Corrected Reads)", paste("Per Litre Filtered"))),
       title="P6") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

#Save all 6 plots in 1
Multiplot_NoZero <- grid.arrange(grobs=list(P1_NoZero,P2_NoZero,P3_NoZero,
                                            P4_NoZero,P5_NoZero,P6_NoZero),
                                 cols=2,
                                 top="Raw Corrected Data")
ggsave(Multiplot_NoZero,
       file="AllSpecies_DataVisualization_Transformations_AllMarkers_NoZero_24Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers
#pull out just the Proportion data for modelling
All_data_NoZero_Prop <- All_data_NoZero_Model %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(All_data_NoZero_Prop$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(All_data_NoZero_Prop$PropCorrectedReadsPerLitre),scientific=F)
#[1] "0.0000007895833" "1.0570512820513"
boxplot(All_data_NoZero_Prop$PropCorrectedReadsPerLitre)
summary(All_data_NoZero_Prop$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.0663942 0.2116392 0.2971252 0.4761385 1.0570513
IQR <- IQR(All_data_NoZero_Prop$PropCorrectedReadsPerLitre)
quartiles <- quantile(All_data_NoZero_Prop$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_Prop_NoOutliers_Temp <- All_data_NoZero_Prop %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(All_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#0 outliers removed

#explanatory variable - Mean Quant score
hist(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
# [1]    0.6029412 248.6284000
summary(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.6029   6.3355  12.1762  18.8796  21.1311 248.6284      117 
boxplot(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(All_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_Prop_NoOutliers <- All_data_NoZero_Prop_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre) 
boxplot(All_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre)
rm(All_data_NoZero_Prop_NoOutliers_Temp)
#182 outliers removed

#still have a proportion over 1, so remove that value
format(range(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F)
All_data_NoZero_Prop_NoOutliers <- All_data_NoZero_Prop_NoOutliers %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 5 rows


#redo plots with no outliers
P1_NoZero_NoOut <- ggplot(All_data_NoZero_Prop_NoOutliers,
                   aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1 No Outliers")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_NoOut <- ggplot(All_data_NoZero_Prop_NoOutliers,
                   aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_NoOut <- ggplot(All_data_NoZero_Prop_NoOutliers,
                   aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc



#Save all 6 plots in 1 (multiplot function from online)
Multiplot_NoZero_NoOut <- grid.arrange(grobs=list(P1_NoZero_NoOut,P1_NoZero,P2_NoZero,
                                                  P3_NoZero_NoOut,P3_NoZero,P4_NoZero,
                                                  P5_NoZero_NoOut,P5_NoZero,P6_NoZero),
                                       cols=3,
                                       top="No Outliers")
ggsave(Multiplot_NoZero_NoOut,
       file="AllSpecies_DataVisualization_Transformations_AllMarkers_NoZero_NoOutliers_24Oct2023.pdf", 
       height=20, width=25,units = "in")


#remove the plots once saved
rm(P1_NoZero_NoOut,P1_NoZero,P2_NoZero,
   P3_NoZero_NoOut,P3_NoZero,P4_NoZero,
   P5_NoZero_NoOut,P5_NoZero,P6_NoZero)

#export data
write.csv(All_data_NoZero_Prop_NoOutliers,
          "AllSpecies_AllMarkers_NoZero_NoOutliers_25Oct2023.csv")


#####All Markers Error Structure Selection #####
#Vars
# River = Random
# Type = Fixed - only 4 types, and should not be nested within river
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed
# Marker = Fixed

#picking error structure
#hist of response variable
hist(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre) 
format(range(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.000001240196" "1.057051282051"
var(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.06758297
mean(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.3001119

descdist(All_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre, boot=500) #likely closest to normal or beta dist

#Normal error distribution
Full_LMM_NoZero_normal_iden <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                       na.action = na.omit,
                                       family =gaussian, #linear
                                       REML=F)
Full_LMM_NoZero_normal_log <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                      na.action = na.omit,
                                      family =gaussian(link="log"), #logarithmic?
                                      REML=F)
Full_LMM_NoZero_normal_logit <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                        na.action = na.omit,
                                        family =gaussian(link="logit"), #logarithmic?
                                        REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the Full model with different error structures, then compare with diagnostic plots

#beta - fall between 0 and 1 (proportion!) 
Full_LMM_NoZero_Beta <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                na.action = na.omit,
                                family = beta_family(),
                                REML=F) #default logit
#beta - fall between 0 and 1 (proportion!) 
Full_LMM_NoZero_BetaBN <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                na.action = na.omit,
                                family = betabinomial(link="logit"),
                                REML=F) #default logit
#Gamma
Full_LMM_NoZero_Gamma <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                 na.action = na.fail,
                                 family=Gamma(link="log"),
                                 REML=F)

#add in negative binomial - can handle zeros?
Full_LMM_NoZero_NB <- glmer.nb(data = All_data_NoZero_Prop_NoOutliers,
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                               na.action = na.omit)

#Model comparisons
#Plots - qqnorm, DHARMa
plot(simulateResiduals(Full_LMM_NoZero_normal_iden))
plot(simulateResiduals(Full_LMM_NoZero_normal_log))
plot(simulateResiduals(Full_LMM_NoZero_normal_logit))
plot(simulateResiduals(Full_LMM_NoZero_Beta))
plot(simulateResiduals(Full_LMM_NoZero_BetaBN)) #not good
plot(simulateResiduals(Full_LMM_NoZero_Gamma))
plot(simulateResiduals(Full_LMM_NoZero_NB)) 

#res vs fit
plot(x=fitted(Full_LMM_NoZero_normal_iden),y=residuals(Full_LMM_NoZero_normal_iden))
plot(x=fitted(Full_LMM_NoZero_normal_log),y=residuals(Full_LMM_NoZero_normal_log))
plot(x=fitted(Full_LMM_NoZero_normal_logit),y=residuals(Full_LMM_NoZero_normal_logit))
plot(x=fitted(Full_LMM_NoZero_Beta),y=residuals(Full_LMM_NoZero_Beta))
plot(x=fitted(Full_LMM_NoZero_Gamma),y=residuals(Full_LMM_NoZero_Gamma))
plot(x=fitted(Full_LMM_NoZero_NB),y=residuals(Full_LMM_NoZero_NB))

#histogram of residuals
hist(residuals(Full_LMM_NoZero_normal_iden),main="Normal")
hist(residuals(Full_LMM_NoZero_normal_log),main="Normal w/log")
hist(residuals(Full_LMM_NoZero_normal_logit),main="Normal w/logit")
hist(residuals(Full_LMM_NoZero_Beta),main="Beta")
hist(residuals(Full_LMM_NoZero_Gamma),main="Gamma")
hist(residuals(Full_LMM_NoZero_NB),main="NB")

#test dispersion
testDispersion(simulateResiduals(Full_LMM_NoZero_normal_iden,1000))
testDispersion(simulateResiduals(Full_LMM_NoZero_normal_log,1000))
testDispersion(simulateResiduals(Full_LMM_NoZero_normal_logit,1000))
testDispersion(simulateResiduals(Full_LMM_NoZero_Beta,1000))
testDispersion(simulateResiduals(Full_LMM_NoZero_Gamma,1000))
testDispersion(simulateResiduals(Full_LMM_NoZero_NB,1000))

#AICc
AICctab(Full_LMM_NoZero_normal_iden,
        Full_LMM_NoZero_normal_log,
        Full_LMM_NoZero_normal_logit,
        Full_LMM_NoZero_Beta,
        Full_LMM_NoZero_BetaBN,
        Full_LMM_NoZero_Gamma,
        Full_LMM_NoZero_NB,
        delta=T,base=T)

#given the data (qq-norm, resid vs fit, AICc), move forward with Beta into variable selection


#####All Markers Variable Selection#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
Full_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                           PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                           na.action = na.omit,
                           family = beta_family(),
                           REML=T)
summary(Full_LMM_NoZero)
AICc(Full_LMM_NoZero) #-943.34

#sort out random using AIC
Random1_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                              PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale,
                              na.action = na.omit,
                              family = beta_family(),
                              REML=T)
Random2_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                              PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run),
                              na.action = na.omit,
                              family = beta_family(),
                              REML=T)
Random3_LMM_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                              PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Code),
                              na.action = na.omit,
                              family = beta_family(),
                              REML=T)

#check AICc of the models
AICc(Full_LMM_NoZero)  #-943.34  
AICc(Random1_LMM_NoZero) #-755.9799
AICc(Random2_LMM_NoZero)  #  -794.9804
AICc(Random3_LMM_NoZero) # -936.4315

#Full Model has the best AIC score

#Select Fixed Effects, use AIC and DO NOT use REML
Full_LMM_NoZero_Fixed <-  glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                                  PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                  na.action = na.fail,
                                  family = beta_family(),
                                  REML=F)
AICc(Full_LMM_NoZero_Fixed) #-974.2258

#use dredge to run all possible sub-models
Dredge_Full_LMM_NoZero_Fixed <- dredge(Full_LMM_NoZero_Fixed)
print(Dredge_Full_LMM_NoZero_Fixed)
# Model selection table 
#     cnd((Int)) dsp((Int)) cnd(DNA) cnd(Mrk) cnd(QMP) cnd(Rsl) cnd(Spc) cnd(Typ) df  logLik   AICc  delta weight
# 30    -0.9274          + -0.25160           0.05270        +        +          11 501.885 -981.4   0.00  0.586
# 32    -0.8181          + -0.25070        +  0.05255        +        +          13 503.180 -979.9   1.54  0.272

#Top 2 models are basically equivalent, based on deltaAIC (<2)
#use the simplest model, with just DNA Concentration,QuantMeanPerLitre,Result,Species


##### AllMarkers Final model#####
FinalModel_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers, 
                             PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|Run) + (1|Code),
                             na.action = na.omit,
                             family = beta_family(),
                             REML=T)
plot(simulateResiduals(FinalModel_NoZero))
summary(FinalModel_NoZero) 
Anova(FinalModel_NoZero,type="III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: PropCorrectedReadsPerLitre
#                       Chisq Df Pr(>Chisq)    
#   (Intercept)       20.351  1  6.446e-06 ***
#   QuantMeanPerLitre 80.576  1  < 2.2e-16 ***
#   DNAConcScale      13.528  1  0.0002351 ***
#   Result            25.703  3  1.101e-05 ***
#   Species           13.197  2  0.0013622 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#plain LM to compare
LM_finalmodel_NoZero <- lm(data = All_data_NoZero_Prop_NoOutliers, 
                           na.action = na.omit,
                           PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + Run + Code)
plot(simulateResiduals(LM_finalmodel_NoZero))
summary(LM_finalmodel_NoZero) 

#Final model with log normal
FinalModel_log_NoZero <- glmmTMB(data = All_data_NoZero_Prop_NoOutliers,
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|Run) + (1|Code),
                                 na.action = na.omit,
                                 REML=T,
                                 family = gaussian(link="logit"))
plot(simulateResiduals(FinalModel_log_NoZero))
summary(FinalModel_log_NoZero) 


AICctab(FinalModel_NoZero,
        LM_finalmodel_NoZero,
        FinalModel_log_NoZero,
        delta=T,base=T)


###plot the model
# Extract the prediction data frame
pred_mm_NoZero <- ggpredict(FinalModel_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model
pred_lm_NoZero <- ggpredict(LM_finalmodel_NoZero, terms = c("QuantMeanPerLitre"))
pred_log_NoZero <- ggpredict(FinalModel_log_NoZero, terms = c("QuantMeanPerLitre"))


# Plot the predictions 
FinalModel_NoZero_plot <- ggplot(pred_mm_NoZero) + 
  geom_line(data = pred_mm_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_mm_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "All Species Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalModel_NoZero_plot

FinalLModel_NoZero_plot <- ggplot(pred_lm_NoZero) + 
  geom_line(data = pred_lm_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_lm_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "LM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + Run + River") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLModel_NoZero_plot

FinalLogModel_NoZero_plot <- ggplot(pred_log_NoZero) + 
  geom_line(data = pred_log_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_log_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Log Normal MM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|Run) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLogModel_NoZero_plot

#export plots
ggsave(FinalModel_NoZero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_NoZero_NoOutliers_GLMMBeta_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLModel_NoZero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_NoZero_NoOutliers_LM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLogModel_NoZero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_NoZero_NoOutliers_LogNormalMM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


###################Atlantic Salmon################
####Data Prep####
#Normalize data by Vol filtered and remove any rows with 0 corrected reads
AtlSalmon_data_NoZero <- AtlSalmon_data_raw %>% 
  filter(CorrectedReads > 0) %>% 
  mutate(RawReadsPerLitre = RawReads/VolFiltered,
         PropRawReadsPerLitre = RawPropReads/VolFiltered,
         CorrectedReadsPerLitre = CorrectedReads/VolFiltered,
         PropCorrectedReadsPerLitre = CorrectedPropReads/VolFiltered,
         QuantMeanPerLitre = QuantMean/VolFiltered) %>% 
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
AtlSalmon_data_NoZero_Model <- AtlSalmon_data_NoZero %>% 
  dplyr::select("SampleID","Type","Name","Code","Date","Marker","Taxon",
                "TotalVertReadsPerSample","RawReads","RawReadsPerLitre",
                "RawPropReads","PropRawReadsPerLitre",
                "CorrectedReads","CorrectedReadsPerLitre",
                "CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFiltered","DNAConc","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result")



#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6


P1_NoZero_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Model,
                              aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +#colour by marker type 
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P2_NoZero_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Model,
                              aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Corrected Reads Per Litre Filtered",
       title="P2")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Model,
                              aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P4_NoZero_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Model,
                              aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Corrected Reads Per Litre Filtered",
       title="P4") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Model,
                              aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P6_NoZero_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Model,
                              aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Corrected Reads)", paste("Per Litre Filtered"))),
       title="P6") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

#Save all 6 plots in 1
Multiplot_NoZero_AtlSalmon <- grid.arrange(grobs=list(P1_NoZero_AtlSalmon,P2_NoZero_AtlSalmon,P3_NoZero_AtlSalmon,
                                                      P4_NoZero_AtlSalmon,P5_NoZero_AtlSalmon,P6_NoZero_AtlSalmon),
                                           cols=2,
                                           top="Raw Corrected Data")
ggsave(Multiplot_NoZero_AtlSalmon,
       file="AtlanticSalmon_DataVisualization_Transformations_AllMarkers_NoZero_24Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers
#pull out just the Proportion data for modelling
AtlSalmon_data_NoZero_Prop <- AtlSalmon_data_NoZero_Model %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre),scientific=F)
#[1] " "0.0000009154" "1.0570512821"
boxplot(AtlSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre)
summary(AtlSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000009 0.0793825 0.2124792 0.2659394 0.4078684 1.0570513 
IQR <- IQR(AtlSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_NoOutliers_Temp <- AtlSalmon_data_NoZero_Prop %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#9 outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
# [1]     1.138026 248.628400
summary(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.138   7.460  12.600  17.438  20.897 248.628      61 
boxplot(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_NoOutliers <- AtlSalmon_data_NoZero_Prop_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre)
rm(AtlSalmon_data_NoZero_Prop_NoOutliers_Temp)
#98 outliers removed


#redo plots with no outliers

P1_NoZero_NoOut_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Prop_NoOutliers,
                                    aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1 No Outliers")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_NoOut_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Prop_NoOutliers,
                                    aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_NoOut_AtlSalmon <- ggplot(AtlSalmon_data_NoZero_Prop_NoOutliers,
                                    aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc



#Save all 6 plots in 1 (multiplot function from online)
Multiplot_NoZero_NoOut_AtlSalmon <- grid.arrange(grobs=list(P1_NoZero_NoOut_AtlSalmon,P1_NoZero_AtlSalmon,P2_NoZero_AtlSalmon,
                                                            P3_NoZero_NoOut_AtlSalmon,P3_NoZero_AtlSalmon,P4_NoZero_AtlSalmon,
                                                            P5_NoZero_NoOut_AtlSalmon,P5_NoZero_AtlSalmon,P6_NoZero_AtlSalmon),
                                                 cols=3,
                                                 top="No Outliers")
ggsave(Multiplot_NoZero_NoOut_AtlSalmon,
       file="AtlanticSalmon_DataVisualization_Transformations_AllMarkers_NoZero_NoOutliers_24Oct2023.pdf", 
       height=20, width=25,units = "in")
#remove the plots once saved
rm(P1_NoZero_NoOut_AtlSalmon,P1_NoZero_AtlSalmon,P2_NoZero_AtlSalmon,
   P3_NoZero_NoOut_AtlSalmon,P3_NoZero_AtlSalmon,P4_NoZero_AtlSalmon,
   P5_NoZero_NoOut_AtlSalmon,P5_NoZero_AtlSalmon,P6_NoZero_AtlSalmon)

#export data
write.csv(AtlSalmon_data_NoZero_Prop_NoOutliers,
          "AtlSalmon_AllMarkers_NoZero_NoOutliers_25Oct2023.csv")

#####All Markers Error Structure Selection #####
#Vars
# River = Random
# Type = Fixed - only 4 types, and should not be nested within river
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed
# Marker = Fixed

#picking error structure
#hist of response variable
hist(AtlSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre) 
format(range(AtlSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.000001240196" "0.864516129032"
var(AtlSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.04284117
mean(AtlSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.2640491

descdist(AtlSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre, boot=500) #likely closest to normal or beta dist

#Normal error distribution
AtlSalmon_Full_LMM_NoZero_normal_iden <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                                       na.action = na.omit,
                                       family =gaussian, #linear
                                       REML=F)
AtlSalmon_Full_LMM_NoZero_normal_log <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                                      na.action = na.omit,
                                      family =gaussian(link="log"), #logarithmic?
                                      REML=F)
AtlSalmon_Full_LMM_NoZero_normal_logit <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                                        na.action = na.omit,
                                        family =gaussian(link="logit"), #logarithmic?
                                        REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the AtlSalmon_Full model with different error structures, then compare with diagnostic plots

#beta - fall between 0 and 1 (proportion!) 
AtlSalmon_Full_LMM_NoZero_Beta <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                                na.action = na.omit,
                                family = beta_family()) #default logit

#Gamma
AtlSalmon_Full_LMM_NoZero_Gamma <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                                 na.action = na.omit,
                                 family=Gamma)

#add in negative binomial - can handle zeros?
AtlSalmon_Full_LMM_NoZero_NB <- glmer.nb(data = AtlSalmon_data_NoZero_Prop_NoOutliers,
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                               na.action = na.omit)


#Model comparisons
#qqnorm DHARMa
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero_normal_iden))
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero_normal_log))
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero_normal_logit))
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero_Beta))
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero_Gamma)) #not great
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero_NB)) 

#res vs fit
plot(x=fitted(AtlSalmon_Full_LMM_NoZero_normal_iden),y=residuals(AtlSalmon_Full_LMM_NoZero_normal_iden))
plot(x=fitted(AtlSalmon_Full_LMM_NoZero_normal_log),y=residuals(AtlSalmon_Full_LMM_NoZero_normal_log))
plot(x=fitted(AtlSalmon_Full_LMM_NoZero_normal_logit),y=residuals(AtlSalmon_Full_LMM_NoZero_normal_logit))
plot(x=fitted(AtlSalmon_Full_LMM_NoZero_Beta),y=residuals(AtlSalmon_Full_LMM_NoZero_Beta))
plot(x=fitted(AtlSalmon_Full_LMM_NoZero_Gamma),y=residuals(AtlSalmon_Full_LMM_NoZero_Gamma))
plot(x=fitted(AtlSalmon_Full_LMM_NoZero_NB),y=residuals(AtlSalmon_Full_LMM_NoZero_NB))


#histogram of residuals
hist(residuals(AtlSalmon_Full_LMM_NoZero_normal_iden),main="Normal")
hist(residuals(AtlSalmon_Full_LMM_NoZero_normal_log),main="Normal w/log")
hist(residuals(AtlSalmon_Full_LMM_NoZero_normal_logit),main="Normal w/logit")
hist(residuals(AtlSalmon_Full_LMM_NoZero_Beta),main="Beta")
hist(residuals(AtlSalmon_Full_LMM_NoZero_Gamma),main="Gamma")
hist(residuals(AtlSalmon_Full_LMM_NoZero_NB),main="NB")

#test dispersion
testDispersion(simulateResiduals(AtlSalmon_Full_LMM_NoZero_normal_iden,1000))
testDispersion(simulateResiduals(AtlSalmon_Full_LMM_NoZero_normal_log,1000))
testDispersion(simulateResiduals(AtlSalmon_Full_LMM_NoZero_normal_logit,1000))
testDispersion(simulateResiduals(AtlSalmon_Full_LMM_NoZero_Beta,1000))
testDispersion(simulateResiduals(AtlSalmon_Full_LMM_NoZero_Gamma,1000))
testDispersion(simulateResiduals(AtlSalmon_Full_LMM_NoZero_NB,1000))

#AICc
AICctab(AtlSalmon_Full_LMM_NoZero_normal_iden,
        AtlSalmon_Full_LMM_NoZero_normal_log,
        AtlSalmon_Full_LMM_NoZero_normal_logit,
        AtlSalmon_Full_LMM_NoZero_Beta,
        AtlSalmon_Full_LMM_NoZero_Gamma,
        AtlSalmon_Full_LMM_NoZero_NB,
        delta=T,base=T)

#given the data, move forward with Beta into variable selection


#####All Markers Variable Selection#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
AtlSalmon_Full_LMM_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                     PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + DNAConcScale + (1|Run) + (1|Code),
                                     na.action = na.omit,
                                     family = beta_family(),
                                     REML=T) 
plot(simulateResiduals(AtlSalmon_Full_LMM_NoZero))
summary(AtlSalmon_Full_LMM_NoZero)
AICc(AtlSalmon_Full_LMM_NoZero) #-860.7958

#sort out random using AIC
AtlSalmon_Random1_LMM_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run),
                                        na.action = na.omit,
                                        family = beta_family(),
                                        REML=T)
AtlSalmon_Random2_LMM_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Code),
                                        na.action = na.omit,
                                        family = beta_family(),
                                        REML=T)
AtlSalmon_Random3_LMM_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result,
                                        na.action = na.omit,
                                        family = beta_family(),
                                        REML=T)

#check AICc of the models
AICc(AtlSalmon_Full_LMM_NoZero)  #-860.7958  
AICc(AtlSalmon_Random1_LMM_NoZero) #-752.7143
AICc(AtlSalmon_Random2_LMM_NoZero)  #  -847.1856
AICc(AtlSalmon_Random3_LMM_NoZero) # -702.9786

#Full Model has the best AIC score

#Select Fixed Effects, use AIC and DO NOT use REML
AtlSalmon_Full_LMM_NoZero_Fixed <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                                           PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                           na.action = na.fail,
                                           family = beta_family(),
                                           REML=F)
AICc(AtlSalmon_Full_LMM_NoZero_Fixed) #-888.1634

#use dredge to run all possible sub-models
Dredge_AtlSalmon_Full_LMM_NoZero_Fixed <- dredge(AtlSalmon_Full_LMM_NoZero_Fixed)
print(Dredge_AtlSalmon_Full_LMM_NoZero_Fixed)
# Model selection table 
#     cnd((Int)) dsp((Int)) cnd(DNA) cnd(Mrk) cnd(QMP) cnd(Rsl) cnd(Typ) df  logLik   AICc delta weight
# 16     -1.754          + -0.11470        +  0.04327        +          11 457.714 -893.0  0.00  0.329
# 15     -1.704          +                 +  0.04026        +          10 456.250 -892.2  0.86  0.215
# 7      -1.819          +                 +  0.04462                    7 452.760 -891.3  1.67  0.143
# 8      -1.869          + -0.09714        +  0.04749                    8 453.771 -891.3  1.70  0.141

#Top 4 models are basically equivalent, based on deltaAIC (<2)
#use the simplest model, with just QuantMeanPerLitre + Marker


##### AllMarkers Final model#####
AtlSalmon_FinalModel_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, 
                             PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + (1|Run) + (1|Code),
                             na.action = na.omit,
                             family = beta_family(),
                             REML=T)
plot(simulateResiduals(AtlSalmon_FinalModel_NoZero))
summary(AtlSalmon_FinalModel_NoZero) 
Anova(AtlSalmon_FinalModel_NoZero,type="III") 

#plain LM to compare
AtlSalmon_LM_finalmodel_NoZero <- lm(data = AtlSalmon_data_NoZero_Prop_NoOutliers, na.action = na.omit,
                                     PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Run + Code)
plot(simulateResiduals(AtlSalmon_LM_finalmodel_NoZero))
summary(AtlSalmon_LM_finalmodel_NoZero) 

#Final model with log normal
AtlSalmon_FinalModel_log_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_NoOutliers, na.action = na.omit,REML=T,family = gaussian(link="logit"),
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + (1|Run) + (1|Code))
plot(simulateResiduals(AtlSalmon_FinalModel_log_NoZero))
summary(AtlSalmon_FinalModel_log_NoZero) 


#AICc
AICctab(AtlSalmon_FinalModel_NoZero,
        AtlSalmon_LM_finalmodel_NoZero,
        AtlSalmon_FinalModel_log_NoZero,delta=T,base=T)

###plot the model
# Extract the prediction data frame
AtlSalmon_pred_mm_NoZero <- ggpredict(AtlSalmon_FinalModel_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model
AtlSalmon_pred_lm_NoZero <- ggpredict(AtlSalmon_LM_finalmodel_NoZero, terms = c("QuantMeanPerLitre"))
AtlSalmon_pred_log_NoZero <- ggpredict(AtlSalmon_FinalModel_log_NoZero, terms = c("QuantMeanPerLitre"))


# Plot the predictions 
AtlSalmon_FinalModel_NoZero_plot <- ggplot(AtlSalmon_pred_mm_NoZero) + 
  geom_line(data = AtlSalmon_pred_mm_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = AtlSalmon_pred_mm_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Atlantic Salmon Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + (1|qPCRRun) + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right",axis.text = element_text(size = 20), axis.title=element_text(size = 20), panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"), panel.grid.major=element_blank());AtlSalmon_FinalModel_NoZero_plot

AtlSalmon_FinalLModel_NoZero_plot <- ggplot(AtlSalmon_pred_lm_NoZero) + 
  geom_line(data = AtlSalmon_pred_lm_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = AtlSalmon_pred_lm_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "LM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + qPCRRun + River") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right",axis.text = element_text(size = 20), axis.title=element_text(size = 20), panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"), panel.grid.major=element_blank());AtlSalmon_FinalLModel_NoZero_plot

AtlSalmon_FinalLogModel_NoZero_plot <- ggplot(AtlSalmon_pred_log_NoZero) + 
  geom_line(data = AtlSalmon_pred_log_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = AtlSalmon_pred_log_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Log Normal MM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + (1|qPCRRun) + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right",axis.text = element_text(size = 20), axis.title=element_text(size = 20), panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"), panel.grid.major=element_blank());AtlSalmon_FinalLogModel_NoZero_plot

#export plots
ggsave(AtlSalmon_FinalModel_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_AllMarkers_NoZero_NoOutliers_GLMMBeta_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(AtlSalmon_FinalLModel_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_AllMarkers_NoZero_NoOutliers_LM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(AtlSalmon_FinalLogModel_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_AllMarkers_NoZero_NoOutliers_LogNormalMM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


###################Arctic Charr################
####Data Prep####
#Normalize data by Vol filtered and remove any rows with 0 corrected reads
ArcCharr_data_NoZero <- ArcticCharr_data_raw %>% 
  filter(CorrectedReads > 0) %>% 
  mutate(RawReadsPerLitre = RawReads/VolFiltered,
         PropRawReadsPerLitre = RawPropReads/VolFiltered,
         CorrectedReadsPerLitre = CorrectedReads/VolFiltered,
         PropCorrectedReadsPerLitre = CorrectedPropReads/VolFiltered,
         QuantMeanPerLitre = QuantMean/VolFiltered) %>% 
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
ArcCharr_data_NoZero_Model <- ArcCharr_data_NoZero %>% 
  dplyr::select("SampleID","Type","Name","Code","Date","Marker","Taxon",
                "TotalVertReadsPerSample","RawReads","RawReadsPerLitre",
                "RawPropReads","PropRawReadsPerLitre",
                "CorrectedReads","CorrectedReadsPerLitre",
                "CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFiltered","DNAConc","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result")



#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6


P1_NoZero_ArcCharr <- ggplot(ArcCharr_data_NoZero_Model,
                             aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))  +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P2_NoZero_ArcCharr <- ggplot(ArcCharr_data_NoZero_Model,
                             aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Corrected Reads Per Litre Filtered",
       title="P2")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_ArcCharr <- ggplot(ArcCharr_data_NoZero_Model,
                             aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P4_NoZero_ArcCharr <- ggplot(ArcCharr_data_NoZero_Model,
                             aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Corrected Reads Per Litre Filtered",
       title="P4") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_ArcCharr <- ggplot(ArcCharr_data_NoZero_Model,
                             aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P6_NoZero_ArcCharr <- ggplot(ArcCharr_data_NoZero_Model,
                             aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Corrected Reads)", paste("Per Litre Filtered"))),
       title="P6") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

#Save all 6 plots in 1
Multiplot_NoZero_ArcCharr <- grid.arrange(grobs=list(P1_NoZero_ArcCharr,P2_NoZero_ArcCharr,P3_NoZero_ArcCharr,
                                                     P4_NoZero_ArcCharr,P5_NoZero_ArcCharr,P6_NoZero_ArcCharr),
                                          cols=2,
                                          top="Raw Corrected Data")
ggsave(Multiplot_NoZero_ArcCharr,
       file="ArcticCharr_DataVisualization_Transformations_AllMarkers_NoZero_24Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers
#pull out just the Proportion data for modelling
ArcCharr_data_NoZero_Prop <- ArcCharr_data_NoZero_Model %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre),scientific=F)
#[1] "0.0000007895833" "1.0101010101010"
boxplot(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre)
summary(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.0376022 0.3781373 0.4396831 0.8547059 1.0101010 
IQR <- IQR(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre)
quartiles <- quantile(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

ArcCharr_data_NoZero_Prop_NoOutliers_Temp <- ArcCharr_data_NoZero_Prop %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper) #remove outliers based on quartiles
rm(IQR, quartiles,Lower,Upper)
hist(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#0 outliers removed

#explanatory variable - Mean Quant score
hist(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
# [1]    0.6029412 170.2106667
summary(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.6029   3.5812   9.8693  27.2356  31.3055 170.2107       47 
boxplot(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(ArcCharr_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

ArcCharr_data_NoZero_Prop_NoOutliers <- ArcCharr_data_NoZero_Prop_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(ArcCharr_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre) 
boxplot(ArcCharr_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre)
rm(ArcCharr_data_NoZero_Prop_NoOutliers_Temp)
#65 outliers removed


#still have a proportion over 1, so remove that value
format(range(ArcCharr_data_NoZero_Prop$PropCorrectedReadsPerLitre),scientific=F)
ArcCharr_data_NoZero_Prop_NoOutliers <- ArcCharr_data_NoZero_Prop_NoOutliers %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 3 rows

#redo plots with no outliers
P1_NoZero_NoOut_ArcCharr <- ggplot(ArcCharr_data_NoZero_Prop_NoOutliers,
                                   aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1 No Outliers")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_NoOut_ArcCharr <- ggplot(ArcCharr_data_NoZero_Prop_NoOutliers,
                                   aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_NoOut_ArcCharr <- ggplot(ArcCharr_data_NoZero_Prop_NoOutliers,
                                   aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc



#Save all 9 plots in 1 (multiplot function from online)
Multiplot_NoZero_NoOut_ArcCharr <- grid.arrange(grobs=list(P1_NoZero_NoOut_ArcCharr,P1_NoZero_ArcCharr,P2_NoZero_ArcCharr,
                                                           P3_NoZero_NoOut_ArcCharr,P3_NoZero_ArcCharr,P4_NoZero_ArcCharr,
                                                           P5_NoZero_NoOut_ArcCharr,P5_NoZero_ArcCharr,P6_NoZero_ArcCharr),
                                                cols=3,
                                                top="No Outliers")
ggsave(Multiplot_NoZero_NoOut_ArcCharr,
       file="ArcticCharr_DataVisualization_Transformations_AllMarkers_NoZero_NoOutliers_24Oct2023.pdf", 
       height=20, width=25,units = "in")

#remove the plots once saved
rm(P1_NoZero_NoOut_ArcCharr,P1_NoZero_ArcCharr,P2_NoZero_ArcCharr,
   P3_NoZero_NoOut_ArcCharr,P3_NoZero_ArcCharr,P4_NoZero_ArcCharr,
   P5_NoZero_NoOut_ArcCharr,P5_NoZero_ArcCharr,P6_NoZero_ArcCharr)

#export data
write.csv(ArcCharr_data_NoZero_Prop_NoOutliers,
          "ArcticCharr_AllMarkers_NoZero_NoOutliers_25Oct2023.csv")


#####All Markers Error Structure Selection #####
#Vars
# River = Random
# Type = Fixed - only 4 types, and should not be nested within river
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed
# Marker = Fixed

#picking error structure
#hist of response variable
hist(ArcCharr_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre) 
format(range(ArcCharr_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.000001359596" "0.999900000000"
var(ArcCharr_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.1338451
mean(ArcCharr_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.5011506

descdist(ArcCharr_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre, boot=500) #likely closest to normal or beta dist

#Normal error distribution
ArcCharr_Full_LMM_NoZero_normal_iden <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                             PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                             na.action = na.omit,
                                             family =gaussian, #linear
                                             REML=F)
ArcCharr_Full_LMM_NoZero_normal_log <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                            na.action = na.omit,
                                            family =gaussian(link="log"), #logarithmic?
                                            REML=F)
ArcCharr_Full_LMM_NoZero_normal_logit <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                              PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                              na.action = na.omit,
                                              family =gaussian(link="logit"), #logarithmic?
                                              REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the Full model with different error structures, then compare with diagnostic plots

#beta - fall between 0 and 1 (proportion!) 
ArcCharr_Full_LMM_NoZero_Beta <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                      na.action = na.omit,
                                      family = beta_family(),
                                      REML=F) #default logit

#Gamma
ArcCharr_Full_LMM_NoZero_Gamma <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                       na.action = na.omit,
                                       family=Gamma)

#add in negative binomial - can handle zeros?
ArcCharr_Full_LMM_NoZero_NB <- glmer.nb(data = ArcCharr_data_NoZero_Prop_NoOutliers,
                                     PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                     na.action = na.omit)

#Model comparisons
#Plots - res vs fit and qqnorm
plot(simulateResiduals(ArcCharr_Full_LMM_NoZero_normal_iden))
plot(simulateResiduals(ArcCharr_Full_LMM_NoZero_normal_log))
plot(simulateResiduals(ArcCharr_Full_LMM_NoZero_normal_logit))
plot(simulateResiduals(ArcCharr_Full_LMM_NoZero_Beta))
plot(simulateResiduals(ArcCharr_Full_LMM_NoZero_Gamma)) #not great
plot(simulateResiduals(ArcCharr_Full_LMM_NoZero_NB)) 

plot(x=fitted(ArcCharr_Full_LMM_NoZero_normal_iden),y=residuals(ArcCharr_Full_LMM_NoZero_normal_iden))
plot(x=fitted(ArcCharr_Full_LMM_NoZero_normal_log),y=residuals(ArcCharr_Full_LMM_NoZero_normal_log))
plot(x=fitted(ArcCharr_Full_LMM_NoZero_normal_logit),y=residuals(ArcCharr_Full_LMM_NoZero_normal_logit))
plot(x=fitted(ArcCharr_Full_LMM_NoZero_Beta),y=residuals(ArcCharr_Full_LMM_NoZero_Beta))
plot(x=fitted(ArcCharr_Full_LMM_NoZero_Gamma),y=residuals(ArcCharr_Full_LMM_NoZero_Gamma))
plot(x=fitted(ArcCharr_Full_LMM_NoZero_NB),y=residuals(ArcCharr_Full_LMM_NoZero_NB))


#histogram of residuals
hist(residuals(ArcCharr_Full_LMM_NoZero_normal_iden),main="Normal")
hist(residuals(ArcCharr_Full_LMM_NoZero_normal_log),main="Normal w/log")
hist(residuals(ArcCharr_Full_LMM_NoZero_normal_logit),main="Normal w/logit")
hist(residuals(ArcCharr_Full_LMM_NoZero_Beta),main="Beta")
hist(residuals(ArcCharr_Full_LMM_NoZero_Gamma),main="Gamma")
hist(residuals(ArcCharr_Full_LMM_NoZero_NB),main="NB")

#test dispersion
testDispersion(simulateResiduals(ArcCharr_Full_LMM_NoZero_normal_iden,1000))
testDispersion(simulateResiduals(ArcCharr_Full_LMM_NoZero_normal_log,1000))
testDispersion(simulateResiduals(ArcCharr_Full_LMM_NoZero_normal_logit,1000))
testDispersion(simulateResiduals(ArcCharr_Full_LMM_NoZero_Beta,1000))
testDispersion(simulateResiduals(ArcCharr_Full_LMM_NoZero_Gamma,1000))
testDispersion(simulateResiduals(ArcCharr_Full_LMM_NoZero_NB,1000))


#AICc
AICctab(ArcCharr_Full_LMM_NoZero_normal_iden,
        ArcCharr_Full_LMM_NoZero_normal_log,
        ArcCharr_Full_LMM_NoZero_normal_logit,
        ArcCharr_Full_LMM_NoZero_Beta,
        ArcCharr_Full_LMM_NoZero_Gamma,
        ArcCharr_Full_LMM_NoZero_NB,
        delta=T,base=T)

#given the data, move forward with Beta into variable selection


#####All Markers Variable Selection#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
ArcCharr_Full_LMM_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                    PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                    na.action = na.omit,
                                    family = beta_family(),
                                    REML=T)
summary(ArcCharr_Full_LMM_NoZero)
AICc(ArcCharr_Full_LMM_NoZero) #-142.1235

#sort out random using AIC
ArcCharr_Random1_LMM_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run),
                                       na.action = na.omit,
                                    family = beta_family(),
                                    REML=T)
ArcCharr_Random2_LMM_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Code),
                                       na.action = na.omit,
                                    family = beta_family(),
                                    REML=T)
ArcCharr_Random3_LMM_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result,
                                       na.action = na.omit,
                                    family = beta_family(),
                                    REML=T)

#check AICc of the models
AICc(ArcCharr_Full_LMM_NoZero)  #-142.1235  
AICc(ArcCharr_Random1_LMM_NoZero) #-121.0967
AICc(ArcCharr_Random2_LMM_NoZero)  #  -144.1206
AICc(ArcCharr_Random3_LMM_NoZero) # -118.1839

#Model2 has the best AIC score - PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Code),

#Select Fixed Effects, use AIC and DO NOT use REML
ArcCharr_Full_LMM_NoZero_Fixed <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Code),
                                       na.action = na.fail,
                                       family = beta_family(),
                                       REML=F)
AICc(ArcCharr_Full_LMM_NoZero_Fixed) #-155.6793

#use dredge to run all possible sub-models
Dredge_ArcCharr_Full_LMM_NoZero_Fixed <- dredge(ArcCharr_Full_LMM_NoZero_Fixed)
print(Dredge_ArcCharr_Full_LMM_NoZero_Fixed)
# Model selection table 
#     cnd((Int)) dsp((Int)) cnd(DNA) cnd(Mrk) cnd(QMP) cnd(Rsl) cnd(Typ) df logLik   AICc delta weight
# 8   -0.439100          +  -0.9122        +  0.02196                    6 87.591 -162.6  0.00  0.273
# 6   -0.610700          +  -0.9411           0.02147                    5 86.077 -161.7  0.85  0.179
# 16  -0.174900          +  -0.8920        +  0.01533        +           8 89.050 -161.0  1.55  0.126
# 12   0.083380          +  -0.9231        +                 +           7 87.710 -160.6  1.98  0.101
#Top 4 models are basically equivalent, based on deltaAIC (deltaAIC<2)
#use the simplest model, with just DNAConcScale + QuantMeanPerLitre


##### AllMarkers Final model#####
ArcCharr_FinalModel_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers,
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|Code),
                                      na.action = na.omit,
                                      family = beta_family(),
                                      REML=T)
plot(simulateResiduals(ArcCharr_FinalModel_NoZero)) 
summary(ArcCharr_FinalModel_NoZero) 
Anova(ArcCharr_FinalModel_NoZero,type="III") 

#plain LM to compare
LM_Charr_finalmodel_NoZero <- lm(data = ArcCharr_data_NoZero_Prop_NoOutliers, na.action = na.omit,
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Code)
plot(simulateResiduals(LM_Charr_finalmodel_NoZero)) 
summary(LM_Charr_finalmodel_NoZero) 

#Final model with log normal
Charr_FinalModel_log_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers,
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|Code),
                                       na.action = na.omit,
                                       REML=T,
                                       family = gaussian(link="logit"))
plot(simulateResiduals(Charr_FinalModel_log_NoZero)) 
summary(Charr_FinalModel_log_NoZero) 



#AICc
AICctab(ArcCharr_FinalModel_NoZero,
        LM_Charr_finalmodel_NoZero,
        Charr_FinalModel_log_NoZero,delta=T,base=T)

###plot the model
# Extract the prediction data frame
pred_mm_Charr_NoZero <- ggpredict(ArcCharr_FinalModel_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model
pred_lm_Charr_NoZero <- ggpredict(LM_Charr_finalmodel_NoZero, terms = c("QuantMeanPerLitre"))
pred_log_Charr_NoZero <- ggpredict(Charr_FinalModel_log_NoZero, terms = c("QuantMeanPerLitre"))


# Plot the predictions 
ArcCharr_FinalModel_NoZero_plot <- ggplot(pred_mm_Charr_NoZero) +
  coord_cartesian(ylim=c(0, 1))+
  geom_line(data = pred_mm_Charr_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_mm_Charr_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = ArcCharr_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Arctic Charr Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())

ArcCharr_FinalLModel_NoZero_plot <- ggplot(pred_lm_Charr_NoZero) +
  coord_cartesian(ylim=c(0, 1))+
  geom_line(data = pred_lm_Charr_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_lm_Charr_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = ArcCharr_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "LM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + qPCRRun + River") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())

ArcCharr_FinalLogModel_NoZero_plot <- ggplot(pred_log_Charr_NoZero) + 
  coord_cartesian(ylim=c(0, 1))+
  geom_line(data = pred_log_Charr_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_log_Charr_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = ArcCharr_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Log Normal MM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + (1|qPCRRun) + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())

#export plots
ggsave(ArcCharr_FinalModel_NoZero_plot, #plot you want to save
       file = "ArcticCharr_AllMarkers_NoZero_GLMMBeta_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(ArcCharr_FinalLModel_NoZero_plot, #plot you want to save
       file = "ArcticCharr_AllMarkers_NoZero_LM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(ArcCharr_FinalLogModel_NoZero_plot, #plot you want to save
       file = "ArcticCharr_AllMarkers_NoZero_LogNormalMM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width






###################Pink Salmon################
####Data Prep####
#Normalize data by Vol filtered and remove any rows with 0 corrected reads
PinkSalmon_data_NoZero <- PinkSalmon_data_raw %>% 
  filter(CorrectedReads > 0) %>% 
  mutate(RawReadsPerLitre = RawReads/VolFiltered,
         PropRawReadsPerLitre = RawPropReads/VolFiltered,
         CorrectedReadsPerLitre = CorrectedReads/VolFiltered,
         PropCorrectedReadsPerLitre = CorrectedPropReads/VolFiltered,
         QuantMeanPerLitre = QuantMean/VolFiltered) %>% 
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
PinkSalmon_data_NoZero_Model <- PinkSalmon_data_NoZero %>% 
  dplyr::select("SampleID","Type","Name","Code","Date","Marker","Taxon",
                "TotalVertReadsPerSample","RawReads","RawReadsPerLitre",
                "RawPropReads","PropRawReadsPerLitre",
                "CorrectedReads","CorrectedReadsPerLitre",
                "CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFiltered","DNAConc","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result")



#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6


P1_NoZero_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Model,
                               aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P2_NoZero_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Model,
                               aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Corrected Reads Per Litre Filtered",
       title="P2")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Model,
                               aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P4_NoZero_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Model,
                               aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Corrected Reads Per Litre Filtered",
       title="P4") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Model,
                               aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc


P6_NoZero_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Model,
                               aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Corrected Reads)", paste("Per Litre Filtered"))),
       title="P6") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

#Save all 6 plots in 1
Multiplot_NoZero_PinkSalmon <- grid.arrange(grobs=list(P1_NoZero_PinkSalmon,P2_NoZero_PinkSalmon,P3_NoZero_PinkSalmon,
                                                       P4_NoZero_PinkSalmon,P5_NoZero_PinkSalmon,P6_NoZero_PinkSalmon),
                                            cols=2,
                                            top="Raw Corrected Data")
ggsave(Multiplot_NoZero_PinkSalmon,
       file="PinkSalmon_DataVisualization_Transformations_AllMarkers_NoZero_25Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers
#pull out just the Proportion data for modelling
PinkSalmon_data_NoZero_Prop <- PinkSalmon_data_NoZero_Model %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(PinkSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(PinkSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre),scientific=F)
#[1] "0.001311111" "0.553541667"
boxplot(PinkSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre)
summary(PinkSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001311 0.009895 0.039588 0.071444 0.094148 0.553542
IQR <- IQR(PinkSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre)
quartiles <- quantile(PinkSalmon_data_NoZero_Prop$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

PinkSalmon_data_NoZero_Prop_NoOutliers_Temp <- PinkSalmon_data_NoZero_Prop %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper) #remove outliers based on quartiles
rm(IQR, quartiles,Lower,Upper)
hist(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#1 outlier removed

#explanatory variable - Mean Quant score
hist(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
# [1]   0.6974047 4.1793307
summary(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.6974  1.5644  1.9305  2.0919  2.3790  4.1793       8 
boxplot(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

PinkSalmon_data_NoZero_Prop_NoOutliers <- PinkSalmon_data_NoZero_Prop_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(PinkSalmon_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre) 
boxplot(PinkSalmon_data_NoZero_Prop_NoOutliers$QuantMeanPerLitre)
rm(PinkSalmon_data_NoZero_Prop_NoOutliers_Temp)
#11 outliers removed

#redo plots with no outliers

P1_NoZero_NoOut_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Prop_NoOutliers,
                                     aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1 No Outliers")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoZero_NoOut_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Prop_NoOutliers,
                                     aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoZero_NoOut_PinkSalmon <- ggplot(PinkSalmon_data_NoZero_Prop_NoOutliers,
                                     aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc



#Save all 6 plots in 1 (multiplot function from online)
Multiplot_NoZero_NoOut_PinkSalmon <- grid.arrange(grobs=list(P1_NoZero_NoOut_PinkSalmon,P1_NoZero_PinkSalmon,P2_NoZero_PinkSalmon,
                                                             P3_NoZero_NoOut_PinkSalmon,P3_NoZero_PinkSalmon,P4_NoZero_PinkSalmon,
                                                             P5_NoZero_NoOut_PinkSalmon,P5_NoZero_PinkSalmon,P6_NoZero_PinkSalmon),
                                                  cols=3,
                                                  top="No Outliers")
ggsave(Multiplot_NoZero_NoOut_PinkSalmon,
       file="PinkSalmon_DataVisualization_Transformations_AllMarkers_NoZero_NoOutliers_25Oct2023.pdf", 
       height=20, width=25,units = "in")


#remove the plots once saved
rm(P1_NoZero_NoOut_PinkSalmon,P1_NoZero_PinkSalmon,P2_NoZero_PinkSalmon,
   P3_NoZero_NoOut_PinkSalmon,P3_NoZero_PinkSalmon,P4_NoZero_PinkSalmon,
   P5_NoZero_NoOut_PinkSalmon,P5_NoZero_PinkSalmon,P6_NoZero_PinkSalmon)

#export data
write.csv(PinkSalmon_data_NoZero_Prop_NoOutliers,
          "PinkSalmon_AllMarkers_NoZero_NoOutliers_25Oct2023.csv")


#####All Markers Error Structure Selection #####
#Vars
# River = Random
# Type = Fixed - only 4 types, and should not be nested within river
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed
# Marker = Fixed

#picking error structure
#hist of response variable
hist(PinkSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre) 
format(range(PinkSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.001311111" "0.174059406"
var(PinkSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.002808512
mean(PinkSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.04359631

descdist(PinkSalmon_data_NoZero_Prop_NoOutliers$PropCorrectedReadsPerLitre, boot=500) #likely closest to normal or beta dist

#Normal error distribution
PinkSalmon_Full_LMM_NoZero_normal_iden <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                                  PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                                  na.action = na.omit,
                                                  family =gaussian, #linear
                                                  REML=F)
PinkSalmon_Full_LMM_NoZero_normal_log <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                                 na.action = na.omit,
                                                 family =gaussian(link="log"), #logarithmic?
                                                 REML=F)
PinkSalmon_Full_LMM_NoZero_normal_logit <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                                   PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                                   na.action = na.omit,
                                                   family =gaussian(link="logit"), #logarithmic?
                                                   REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the Full model with different error structures, then compare with diagnostic plots

#beta - fall between 0 and 1 (proportion!) 
PinkSalmon_Full_LMM_NoZero_Beta <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                           PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                           na.action = na.omit,
                                           family = beta_family(),
                                           REML=F) #default logit

#Gamma
PinkSalmon_Full_LMM_NoZero_Gamma <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                            na.action = na.omit,
                                            family=Gamma,
                                            REML=F)

#add in negative binomial - can handle zeros?
PinkSalmon_Full_LMM_NoZero_NB <- glmer.nb(data = PinkSalmon_data_NoZero_Prop_NoOutliers,
                                          PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                          na.action = na.omit)


#Model comparisons
#Plots - qqnorm DHARMa
plot(simulateResiduals(PinkSalmon_Full_LMM_NoZero_normal_iden))
plot(simulateResiduals(PinkSalmon_Full_LMM_NoZero_normal_log))
plot(simulateResiduals(PinkSalmon_Full_LMM_NoZero_normal_logit))
plot(simulateResiduals(PinkSalmon_Full_LMM_NoZero_Beta))
plot(simulateResiduals(PinkSalmon_Full_LMM_NoZero_Gamma))
plot(simulateResiduals(PinkSalmon_Full_LMM_NoZero_NB)) 

#res vs fit
plot(x=fitted(PinkSalmon_Full_LMM_NoZero_normal_iden),y=residuals(PinkSalmon_Full_LMM_NoZero_normal_iden))
plot(x=fitted(PinkSalmon_Full_LMM_NoZero_normal_log),y=residuals(PinkSalmon_Full_LMM_NoZero_normal_log))
plot(x=fitted(PinkSalmon_Full_LMM_NoZero_normal_logit),y=residuals(PinkSalmon_Full_LMM_NoZero_normal_logit))
plot(x=fitted(PinkSalmon_Full_LMM_NoZero_Beta),y=residuals(PinkSalmon_Full_LMM_NoZero_Beta))
plot(x=fitted(PinkSalmon_Full_LMM_NoZero_Gamma),y=residuals(PinkSalmon_Full_LMM_NoZero_Gamma))
plot(x=fitted(PinkSalmon_Full_LMM_NoZero_NB),y=residuals(PinkSalmon_Full_LMM_NoZero_NB))

#histogram of residuals
hist(residuals(PinkSalmon_Full_LMM_NoZero_normal_iden),main="Normal")
hist(residuals(PinkSalmon_Full_LMM_NoZero_normal_log),main="Normal w/log")
hist(residuals(PinkSalmon_Full_LMM_NoZero_normal_logit),main="Normal w/logit")
hist(residuals(PinkSalmon_Full_LMM_NoZero_Beta),main="Beta")
hist(residuals(PinkSalmon_Full_LMM_NoZero_Gamma),main="Gamma")
hist(residuals(PinkSalmon_Full_LMM_NoZero_NB),main="NB")

#test dispersion
testDispersion(simulateResiduals(PinkSalmon_Full_LMM_NoZero_normal_iden,1000))
testDispersion(simulateResiduals(PinkSalmon_Full_LMM_NoZero_normal_log,1000))
testDispersion(simulateResiduals(PinkSalmon_Full_LMM_NoZero_normal_logit,1000))
testDispersion(simulateResiduals(PinkSalmon_Full_LMM_NoZero_Beta,1000))
testDispersion(simulateResiduals(PinkSalmon_Full_LMM_NoZero_Gamma,1000))
testDispersion(simulateResiduals(PinkSalmon_Full_LMM_NoZero_NB,1000))

#AICc
AICctab(PinkSalmon_Full_LMM_NoZero_normal_iden,
        PinkSalmon_Full_LMM_NoZero_normal_log,
        PinkSalmon_Full_LMM_NoZero_normal_logit,
        PinkSalmon_Full_LMM_NoZero_Beta,
        PinkSalmon_Full_LMM_NoZero_Gamma,
        PinkSalmon_Full_LMM_NoZero_NB,
        delta=T,base=T)


#logit, Beta, and Gamma all have the best AICc
#based on residual vs fit and qqnorm, move forward with beta


#####All Markers Variable Selection#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
PinkSalmon_Full_LMM_NoZero <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run) + (1|Code),
                                      na.action = na.omit,
                                      family = beta_family(),
                                      REML=T) #default logit
summary(PinkSalmon_Full_LMM_NoZero)
AICc(PinkSalmon_Full_LMM_NoZero) #22.41733

#sort out random using AIC
PinkSalmon_Random1_LMM_NoZero <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                         PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Run),
                                         na.action = na.omit,
                                         family = beta_family(),
                                         REML=T)
PinkSalmon_Random2_LMM_NoZero <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                         PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result + (1|Code),
                                         na.action = na.omit,
                                         family = beta_family(),
                                         REML=T)
PinkSalmon_Random3_LMM_NoZero <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                         PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result,
                                         na.action = na.omit,
                                         family = beta_family(),
                                         REML=T)

#AICc
AICctab(PinkSalmon_Full_LMM_NoZero,
        PinkSalmon_Random1_LMM_NoZero,
        PinkSalmon_Random2_LMM_NoZero,
        PinkSalmon_Random3_LMM_NoZero,delta=T,base=T)

#Model3 has the best AIC score - PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result,

#Select Fixed Effects, use AIC and DO NOT use REML
PinkSalmon_Full_LMM_NoZero_Fixed <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers, 
                                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + DNAConcScale + Type + Result,
                                            na.action = na.fail,
                                            family = beta_family(),
                                            REML=F)
AICc(PinkSalmon_Full_LMM_NoZero_Fixed) #-24.65612

#use dredge to run all possible sub-models
Dredge_PinkSalmon_Full_LMM_NoZero_Fixed <- dredge(PinkSalmon_Full_LMM_NoZero_Fixed)
print(Dredge_PinkSalmon_Full_LMM_NoZero_Fixed)
# Model selection table 
#     cnd((Int)) dsp((Int))  cnd(DNA) cnd(Mrk) cnd(QMP) cnd(Rsl) cnd(Typ) df logLik  AICc delta weight
# 9      -3.973          +                                    +           3 34.205 -60.2  0.00  0.263
# 1      -3.086          +                                                2 32.438 -59.9  0.35  0.221
# 5      -4.133          +                      0.5965                    3 33.492 -58.8  1.43  0.129
#Top 3 models are basically equivalent, based on deltaAIC (deltaAIC<2)
#use the simplest model, with just QuantMeanPerLitre


##### AllMarkers Final model#####
PinkSalmon_FinalModel_NoZero <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers,
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre,
                                        na.action = na.omit,
                                        family = beta_family(),
                                        REML=T)
summary(PinkSalmon_FinalModel_NoZero) 
Anova(PinkSalmon_FinalModel_NoZero,type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: PropCorrectedReadsPerLitre
#                       Chisq Df Pr(>Chisq)    
# (Intercept)       24.8203  1  6.293e-07 ***
# QuantMeanPerLitre  1.9418  1     0.1635    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#plain LM to compare
LM_PinkSalmon_finalmodel_NoZero <- lm(data = PinkSalmon_data_NoZero_Prop_NoOutliers, na.action = na.omit,
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre)
summary(LM_PinkSalmon_finalmodel_NoZero) 
anova(LM_PinkSalmon_finalmodel_NoZero)

#Final model with log normal
PinkSalmon_FinalModel_log_NoZero <- glmmTMB(data = PinkSalmon_data_NoZero_Prop_NoOutliers,
                                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre,
                                            na.action = na.omit,
                                            REML=T,
                                            family = gaussian(link="logit"))
summary(PinkSalmon_FinalModel_log_NoZero) 

#AICc
AICctab(PinkSalmon_FinalModel_NoZero,
       LM_PinkSalmon_finalmodel_NoZero,
       PinkSalmon_FinalModel_log_NoZero,
       delta=T,base=T)

###plot the model
# Extract the prediction data frame
pred_mm_PinkSalmon_NoZero <- ggpredict(PinkSalmon_FinalModel_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model
pred_lm_PinkSalmon_NoZero <- ggpredict(LM_PinkSalmon_finalmodel_NoZero, terms = c("QuantMeanPerLitre"))
pred_log_PinkSalmon_NoZero <- ggpredict(PinkSalmon_FinalModel_log_NoZero, terms = c("QuantMeanPerLitre"))


# Plot the predictions 
PinkSalmon_FinalModel_NoZero_plot <- ggplot(pred_mm_PinkSalmon_NoZero) +
  geom_line(data = pred_mm_PinkSalmon_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_mm_PinkSalmon_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = PinkSalmon_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Pink Salmon Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre") +
  theme_bw()+
  theme(legend.position = "right", axis.text = element_text(size = 20),axis.title=element_text(size = 20),panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),panel.grid.major=element_blank());PinkSalmon_FinalModel_NoZero_plot

PinkSalmon_FinalLModel_NoZero_plot <- ggplot(pred_lm_PinkSalmon_NoZero) +
  # coord_cartesian(ylim=c(0,0.2))+
  geom_line(data = pred_lm_PinkSalmon_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_lm_PinkSalmon_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = PinkSalmon_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "LM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre") +
  theme_bw()+
  theme(legend.position = "right", axis.text = element_text(size = 20),axis.title=element_text(size = 20),panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),panel.grid.major=element_blank());PinkSalmon_FinalLModel_NoZero_plot

PinkSalmon_FinalLogModel_NoZero_plot <- ggplot(pred_log_PinkSalmon_NoZero) + 
  geom_line(data = pred_log_PinkSalmon_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_log_PinkSalmon_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = PinkSalmon_data_NoZero_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Log Normal MM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre") +
  theme_bw()+
  theme(legend.position = "right", axis.text = element_text(size = 20),axis.title=element_text(size = 20),panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),panel.grid.major=element_blank());PinkSalmon_FinalLogModel_NoZero_plot

#export plots
ggsave(PinkSalmon_FinalModel_NoZero_plot, #plot you want to save
       file = "PinkSalmon_AllMarkers_NoZero_NoOutliers_GLMMBeta_25Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(PinkSalmon_FinalLModel_NoZero_plot, #plot you want to save
       file = "PinkSalmon_AllMarkers_NoZero_NoOutliers_LM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(PinkSalmon_FinalLogModel_NoZero_plot, #plot you want to save
       file = "PinkSalmon_AllMarkers_NoZero_NoOutliers_LogNormalMM_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width







######Multiplot All FinalModel Plots#####
#Save all 6 plots in 1
Final_Multiplot_NoZero <- grid.arrange(grobs=list(FinalModel_NoZero_plot,
                                                  AtlSalmon_FinalModel_NoZero_plot,
                                                  ArcCharr_FinalModel_NoZero_plot,
                                                  PinkSalmon_FinalModel_NoZero_plot),
                                       cols=2,
                                       top="Beta GLMM")
ggsave(Final_Multiplot_NoZero,
       file="AllSpecies_AllMarkers_NoZero_NoOutliers_GLMMMultiplot_25Oct2023.pdf", 
       height=20, width=25,units = "in")


###############Models with 0 read count included###########
##############All Species################
####Data Prep####
#Normalize data by Vol filtered
All_data <- All_data_raw %>% 
  mutate(RawReadsPerLitre = RawReads/VolFiltered,
         PropRawReadsPerLitre = RawPropReads/VolFiltered,
         CorrectedReadsPerLitre = CorrectedReads/VolFiltered,
         PropCorrectedReadsPerLitre = CorrectedPropReads/VolFiltered,
         QuantMeanPerLitre = QuantMean/VolFiltered) %>% 
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
All_data_Model <- All_data %>% 
  dplyr::select("SampleID","Type","Name","Code","Date","Marker","Taxon",
                "TotalVertReadsPerSample","RawReads","RawReadsPerLitre",
                "RawPropReads","PropRawReadsPerLitre",
                "CorrectedReads","CorrectedReadsPerLitre",
                "CorrectedPropReads","PropCorrectedReadsPerLitre",
                "VolFiltered","DNAConc","DNAConcScale","Run","QuantMean","QuantMeanPerLitre","Result","Species") %>% 
  filter(!is.na(QuantMeanPerLitre)) #remove any rows with an NA in QuantMeanPerLitre

#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6


P1 <- ggplot(All_data_Model,
             aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P2 <- ggplot(All_data_Model,
                    aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Corrected Reads Per Litre Filtered",
       title="P2")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3 <- ggplot(All_data_Model,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())#changes the axes, etc

P4 <- ggplot(All_data_Model,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Corrected Reads Per Litre Filtered",
       title="P4") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5 <- ggplot(All_data_Model,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P6 <- ggplot(All_data_Model,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Corrected Reads)", paste("Per Litre Filtered"))),
       title="P6") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

#Save all 6 plots in 1
Multiplot <- grid.arrange(grobs=list(P1,P2,P3,
                                     P4,P5,P6),
                                 cols=2,
                                 top="Raw Corrected Data")
ggsave(Multiplot,
       file="AllSpecies_DataVisualization_Transformations_AllMarkers_25Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers - NOTE: Data is so skewed removing outliers pulls out most if the data - not good
#pull out just the Proportion data for modelling
All_data_Prop <- All_data_Model %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))


#Response variable - proportion of total reads
hist(All_data_Prop$PropCorrectedReadsPerLitre,breaks = 50) #data is not normal, can't use z-score
format(range(All_data_Prop$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.000000" "1.057051"
boxplot(All_data_Prop$PropCorrectedReadsPerLitre)
summary(All_data_Prop$PropCorrectedReadsPerLitre)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.1443  0.2452  0.4010  1.0571
IQR <- IQR(All_data_Prop$PropCorrectedReadsPerLitre)
quartiles <- quantile(All_data_Prop$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_Prop_NoOutliers_Temp <- All_data_Prop %>%
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre,breaks=50)
boxplot(All_data_Prop_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#5 outliers removed

#explanatory variable - Mean Quant score - Data is so skewed removing outliers removes MOST of the data..don't remove outliers from MeanQuantScore
hist(All_data_Prop_NoOutliers_Temp$QuantMeanPerLitre,breaks=50) #data is not normal, can't use z-score
range(All_data_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
# [1]   0.2315217 345.7539216
summary(All_data_Prop_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.2315   3.9384   9.6377  16.6096  18.0002 345.7539
boxplot(All_data_Prop_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(All_data_Prop_NoOutliers_Temp$QuantMeanPerLitre, na.rm=T)
quartiles <- quantile(All_data_Prop_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_Prop_NoOutliers <- All_data_Prop_NoOutliers_Temp %>%
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_Prop_NoOutliers$QuantMeanPerLitre,breaks=50)
boxplot(All_data_Prop_NoOutliers$QuantMeanPerLitre)
rm(All_data_Prop_NoOutliers_Temp)
# 71 outliers removed

#remove proportions greater than 1
format(range(All_data_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F)
All_data_Prop_NoOutliers <- All_data_Prop_NoOutliers %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 6 rows


#redo plots with no outliers
P1_NoOutliers <- ggplot(All_data_Prop_NoOutliers,
                          aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="P1 No Outliers")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P3_NoOutliers <- ggplot(All_data_Prop_NoOutliers,
                          aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop. Total Reads Per Litre Filtered",
       title="P3 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc

P5_NoOutliers <- ggplot(All_data_Prop_NoOutliers,
                          aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker)) +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop. Total Reads)", paste("Per Litre Filtered"))),
       title="P5 No Outliers") +
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the axes, etc



#Save all 6 plots in 1 (multiplot function from online)
Multiplot_NoOutliers <- grid.arrange(grobs=list(P1_NoOutliers,P1,P2,
                                           P3_NoOutliers,P3,P4,
                                           P5_NoOutliers,P5,P6),
                                       cols=3)
ggsave(Multiplot_NoOutliers,
       file="AllSpecies_DataVisualization_Transformations_AllMarkers_NoOutliers_25Oct2023.pdf", 
       height=20, width=25,units = "in")


#remove the plots once saved
rm(P1_NoOutliers,P1,P2,
   P3_NoOutliers,P3,P4,
   P5_NoOutliers,P5,P6)

#export data
write.csv(All_data_Prop_NoOutliers,
          "AllSpecies_AllMarkers_NoOutliers_25Oct2023.csv")

#####All Markers Error Structure Selection #####
#Vars
# River = Random
# Type = Fixed - only 4 types, and should not be nested within river
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed
# Marker = Fixed

All_data_Prop_NoOutliers

#picking error structure
#hist of response variable
hist(All_data_Prop_NoOutliers$PropCorrectedReadsPerLitre,breaks=50) 
format(range(All_data_Prop_NoOutliers$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.0000" "0.9999"
var(All_data_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.06722745
mean(All_data_Prop_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.2207186

descdist(All_data_Prop_NoOutliers$PropCorrectedReadsPerLitre, boot=500) #likely closest to normal or beta dist

#Normal error distribution
Full_Zero_LMM_normal_iden <- glmmTMB(data = All_data_Prop_NoOutliers, 
                                       PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                       na.action = na.omit,
                                       family =gaussian, #linear
                                       REML=F)
Full_Zero_LMM_normal_log <- glmmTMB(data = All_data_Prop_NoOutliers, 
                                      PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                      na.action = na.omit,
                                      family =gaussian(link="log"), #logarithmic?
                                      REML=F, start=0)
Full_Zero_LMM_normal_logit <- glmmTMB(data = All_data_Prop_NoOutliers, 
                                        PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                        na.action = na.omit,
                                        family =gaussian(link="logit"), #logarithmic?
                                        REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the Full model with different error structures, then compare with diagnostic plots

#Beta distribution was the most appropriate for the data with the zeros removed, but cannot handle zeros in the data
#a zero inflation model can be run, but what defines the 0?

#beta - fall between 0 and 1 (proportion!) 
Full_Zero_LMM_Beta1 <- glmmTMB(data = All_data_Prop_NoOutliers, 
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                               na.action = na.omit,
                               family = beta_family(),#default logit
                               REML=F,
                               ziformula = ~1)  #ziformula=~1 means probability of structural zero is the same for all obs
Full_Zero_LMM_Beta2 <- glmmTMB(data = All_data_Prop_NoOutliers, 
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                               na.action = na.omit,
                               family = beta_family(),#default logit
                               REML=F,
                               ziformula = ~Result) # probability of structural zero varies by Result
Full_Zero_LMM_Beta3 <- glmmTMB(data = All_data_Prop_NoOutliers, 
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                               na.action = na.omit,
                               family = beta_family(),#default logit
                               REML=F,
                               ziformula = ~DNAConcScale) # probability of structural zero varies by DNAConcScale
Full_Zero_LMM_Beta4 <- glmmTMB(data = All_data_Prop_NoOutliers, 
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                               na.action = na.omit,
                               family = beta_family(),#default logit
                               REML=F,
                               ziformula = ~.) # probability of structural zero is the same as the main formula

#hurdle NB from glmmTMB
Full_Zero_LMM_NB <- glmmTMB(data = All_data_Prop_NoOutliers, 
                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                            na.action = na.omit,
                            family = truncated_nbinom1(link="log"),
                            REML=F,
                            ziformula = ~.) 


#Model comparisons
#Plots - qqnorm, DHARMa
plot(simulateResiduals(Full_Zero_LMM_normal_iden)) #bad
plot(simulateResiduals(Full_Zero_LMM_normal_log)) #bad
plot(simulateResiduals(Full_Zero_LMM_normal_logit)) #bad
plot(simulateResiduals(Full_Zero_LMM_Beta1),main="beta1") #bad
plot(simulateResiduals(Full_Zero_LMM_Beta2),main="beta2") #bad
plot(simulateResiduals(Full_Zero_LMM_Beta3),main="beta3") #bad
plot(simulateResiduals(Full_Zero_LMM_Beta4),main="beta4") #pretty good!
plot(simulateResiduals(Full_Zero_LMM_NB))  #very bad

#res vs fit
plot(x=fitted(Full_Zero_LMM_normal_iden),y=residuals(Full_Zero_LMM_normal_iden))
plot(x=fitted(Full_Zero_LMM_normal_log),y=residuals(Full_Zero_LMM_normal_log))
plot(x=fitted(Full_Zero_LMM_normal_logit),y=residuals(Full_Zero_LMM_normal_logit))
plot(x=fitted(Full_Zero_LMM_Beta1),y=residuals(Full_Zero_LMM_Beta1))
plot(x=fitted(Full_Zero_LMM_Beta2),y=residuals(Full_Zero_LMM_Beta2))
plot(x=fitted(Full_Zero_LMM_Beta3),y=residuals(Full_Zero_LMM_Beta3))
plot(x=fitted(Full_Zero_LMM_Beta4),y=residuals(Full_Zero_LMM_Beta4))
plot(x=fitted(Full_Zero_LMM_NB),y=residuals(Full_Zero_LMM_NB))

#histogram of residuals
hist(residuals(Full_Zero_LMM_normal_iden),main="Normal")
hist(residuals(Full_Zero_LMM_normal_log),main="Normal w/log")
hist(residuals(Full_Zero_LMM_normal_logit),main="Normal w/logit")
hist(residuals(Full_Zero_LMM_Beta1),main="Beta1")
hist(residuals(Full_Zero_LMM_Beta2),main="Beta2")
hist(residuals(Full_Zero_LMM_Beta3),main="Beta3")
hist(residuals(Full_Zero_LMM_Beta4),main="Beta4")
hist(residuals(Full_Zero_LMM_NB),main="NB")

#test dispersion
testDispersion(simulateResiduals(Full_Zero_LMM_normal_iden,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_normal_log,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_normal_logit,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_Beta1,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_Beta2,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_Beta3,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_Beta4,1000))
testDispersion(simulateResiduals(Full_Zero_LMM_NB,1000))

#AIC
AICctab(Full_Zero_LMM_normal_iden,
       Full_Zero_LMM_normal_log,
       Full_Zero_LMM_normal_logit,
       Full_Zero_LMM_Beta1,
       Full_Zero_LMM_Beta2,
       Full_Zero_LMM_Beta3,
       Full_Zero_LMM_Beta4,
       Full_Zero_LMM_NB,
       delta=T,base=T)


#given the data (qq-norm, resid vs fit, AICc), move forward with Beta4 into variable selection - the normal dist, while having good AIC, have very bad qqnorm, etc


#####All Markers Variable Selection#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
Full_LMM_Zero <- glmmTMB(data = All_data_Prop_NoOutliers, 
                         PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                         na.action = na.omit,
                         family = beta_family(),#default logit
                         REML=T,
                         ziformula = ~.) 
summary(Full_LMM_Zero)

#sort out random using AIC
Random1_LMM_Zero <- glmmTMB(data = All_data_Prop_NoOutliers, 
                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run),
                            na.action = na.omit,
                            family = beta_family(),#default logit
                            REML=T,
                            ziformula = ~.)
Random2_LMM_Zero <- glmmTMB(data = All_data_Prop_NoOutliers, 
                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Code),
                            na.action = na.omit,
                            family = beta_family(),#default logit
                            REML=T,
                            ziformula = ~.)
Random3_LMM_Zero <- glmmTMB(data = All_data_Prop_NoOutliers, 
                            PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale,
                            na.action = na.omit,
                            family = beta_family(),#default logit
                            REML=T,
                            ziformula = ~.)

#AICc
AICctab(Full_LMM_Zero,
       Random1_LMM_Zero,
       Random2_LMM_Zero,
       Random3_LMM_Zero,
       delta=T, base=T) 

#Full Model has the best AIC score

#Select Fixed Effects, use AIC and DO NOT use REML
Full_LMM_Zero_Fixed <-  glmmTMB(data = All_data_Prop_NoOutliers, 
                                PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + Marker + Type + Result + Species + DNAConcScale + (1|Run) + (1|Code),
                                na.action = na.fail,
                                family = beta_family(),#default logit
                                REML=F,
                                ziformula = ~.)
AIC(Full_LMM_Zero_Fixed) #-22.07563

#use dredge to run all possible sub-models
Dredge_Full_LMM_Zero_Fixed <- dredge(Full_LMM_Zero_Fixed)
print(Dredge_Full_LMM_Zero_Fixed)
# Model selection table 
#       cnd((Int)) zi((Int)) dsp((Int)) cnd(DNA) cnd(Mrk) cnd(QMP) cnd(Rsl) cnd(Spc) cnd(Typ)   zi(DNA) zi(Mrk) zi(QMP) zi(Rsl) zi(Spc) zi(Typ) df   logLik  AICc  delta weight
# 1950    -0.9217  -2.24700          + -0.24610           0.05270        +        +                          + -0.1143       +       +         22   39.585 -34.2   0.00  0.303
# 2014    -0.9217  -2.19000          + -0.24610           0.05270        +        +           0.130700       + -0.1183       +       +         23   40.220 -33.4   0.82  0.201
# 1952    -0.8125  -2.24700          + -0.24510        +  0.05255        +        +                          + -0.1143       +       +         24   40.880 -32.6   1.59  0.137

#Top 3 models are basically equivalent, based on deltaAIC (<2)
#Main: QuantMeanPerLitre + DNAConcScale + Result + Species + (1|Run) + (1|Code)
#ZI: Marker + QuantMeanPerLitre + Result + Species

##### AllMarkers Final model#####
FinalModel_Zero <- glmmTMB(data = All_data_Prop_NoOutliers, 
                           PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|Run) + (1|Code),
                           na.action = na.fail,
                           family = beta_family(),#default logit
                           REML=T,
                           ziformula = ~ QuantMeanPerLitre + Marker + Result + Species)
plot(simulateResiduals(FinalModel_Zero))
summary(FinalModel_Zero) 
Anova(FinalModel_Zero,type="III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: PropCorrectedReadsPerLitre
#                       Chisq Df Pr(>Chisq)    
#   (Intercept)       20.351  1  6.446e-06 ***
#   QuantMeanPerLitre 80.576  1  < 2.2e-16 ***
#   DNAConcScale      13.528  1  0.0002351 ***
#   Result            25.703  3  1.101e-05 ***
#   Species           13.197  2  0.0013622 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#plain LM to compare
LM_finalmodel_Zero <- lm(data = All_data_Prop_NoOutliers, 
                         na.action = na.omit,
                         PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + Run + Code)
plot(simulateResiduals(LM_finalmodel_Zero))
summary(LM_finalmodel_Zero) 

#Final model with log normal
FinalModel_log_Zero <- glmmTMB(data = All_data_Prop_NoOutliers,
                               PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|Run) + (1|Code),
                               na.action = na.omit,
                               REML=T,
                               family = gaussian(link="logit"))
plot(simulateResiduals(FinalModel_log_Zero))
summary(FinalModel_log_Zero) 


AICctab(FinalModel_Zero,
        LM_finalmodel_Zero,
        FinalModel_log_Zero,
        delta=T,base=T)

###plot the model
# Extract the prediction data frame
pred_mm_Zero <- ggpredict(FinalModel_Zero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model
pred_lm_Zero <- ggpredict(LM_finalmodel_Zero, terms = c("QuantMeanPerLitre"))
pred_log_Zero <- ggpredict(FinalModel_log_Zero, terms = c("QuantMeanPerLitre"))


# Plot the predictions 
FinalModel_Zero_plot <- ggplot(pred_mm_Zero) + 
  geom_line(data = pred_mm_Zero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_mm_Zero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "All Species Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalModel_Zero_plot

FinalLModel_Zero_plot <- ggplot(pred_lm_Zero) + 
  geom_line(data = pred_lm_Zero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_lm_Zero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "LM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + Run + River") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLModel_Zero_plot

FinalLogModel_Zero_plot <- ggplot(pred_log_Zero) + 
  geom_line(data = pred_log_Zero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred_log_Zero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = All_data_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Log Normal MM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + Result + Species + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLogModel_Zero_plot

#export plots
ggsave(FinalModel_Zero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_Zero_NoOutliers_GLMMBeta_25Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLModel_Zero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_Zero_NoOutliers_LM_25Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLogModel_Zero_plot, #plot you want to save
       file = "AllSpecies_AllMarkers_Zero_NoOutliers_LogNormalMM_25Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


######save workspace######
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
