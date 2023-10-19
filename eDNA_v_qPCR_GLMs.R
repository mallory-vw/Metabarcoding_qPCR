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
library("gridExtra")
library("car")
library("MuMIn")

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load data (see eDNA_v_qPCR_DataProcessing.R)####
AtlSalmon_data_raw <- read.csv("qPCR_eDNA_SalmoSalar_CleanedData_11Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

PinkSalmon_data_raw <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_11Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

ArcticCharr_data_raw <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_11Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X)) 

####Data Prep####
#calculate QuantMean
AtlSalmon_data_int <- AtlSalmon_data_raw %>% 
  mutate(QuantMean = rowMeans(dplyr::select(.,Quant1,Quant2,Quant3), na.rm=TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)


#Normalize CorrectedReads,CorrectedPercentReads,QuantMean by Vol filtered
AtlSalmon_data <- AtlSalmon_data_int %>% 
  mutate(NormPercReads_12Steleo = CorrectedPercentReads_12Steleo/VolFiltered,
         NormCorrectedReads_12Steleo = CorrectedReads_12Steleo/VolFiltered,
         NormPercReads_FISHE = CorrectedPercentReads_FISHE/VolFiltered,
         NormCorrectedReads_FISHE = CorrectedReads_FISHE/VolFiltered,
         NormPercReads_MIFISHU = CorrectedPercentReads_MIFISHU/VolFiltered,
         NormCorrectedReads_MIFISHU = CorrectedReads_MIFISHU/VolFiltered,
         NormQuantMean = QuantMean/VolFiltered) %>% 
  relocate(NormPercReads_12Steleo, .after=CorrectedPercentReads_12Steleo) %>% 
  relocate(NormCorrectedReads_12Steleo, .after = CorrectedReads_12Steleo) %>% 
  relocate(NormPercReads_FISHE, .after=CorrectedPercentReads_FISHE) %>% 
  relocate(NormCorrectedReads_FISHE, .after = CorrectedReads_FISHE) %>% 
  relocate(NormPercReads_MIFISHU, .after=CorrectedPercentReads_MIFISHU) %>% 
  relocate(NormCorrectedReads_MIFISHU, .after = CorrectedReads_MIFISHU) %>% 
  relocate(NormQuantMean, .after=QuantMean) %>% 
  mutate(RiverCode=factor(RiverCode),
         Type=factor(Type),
         Run=factor(Run),
         Year=factor(Year),
         DNAConcScale=scale(DNAConc_pg_uL, center=T, scale=T)[,1],
         log10DNAConc=log10(DNAConc_pg_uL)) %>% 
  mutate(NormPropReads_12Steleo=NormPercReads_12Steleo/100,
         NormPropReads_FISHE=NormPercReads_FISHE/100,
         NormPropReads_MIFISHU=NormPercReads_MIFISHU/100) %>% 
  relocate(NormPropReads_12Steleo, .after = NormPercReads_12Steleo) %>% 
  relocate(NormPropReads_FISHE, .after = NormPercReads_FISHE) %>% 
  relocate(NormPropReads_MIFISHU, .after = NormPercReads_MIFISHU)

#pivot data into long format
AtlSalmon_data_long <- AtlSalmon_data %>% 
  dplyr::select(SampleID,Type,RiverCode,Year,VolFiltered,
               DNAConc_pg_uL,DNAConcScale,log10DNAConc,
               Run,QuantMean,NormQuantMean,
               CorrectedReads_12Steleo,NormCorrectedReads_12Steleo,CorrectedPercentReads_12Steleo,NormPercReads_12Steleo,NormPropReads_12Steleo,
               CorrectedReads_FISHE,NormCorrectedReads_FISHE,CorrectedPercentReads_FISHE,NormPercReads_FISHE,NormPropReads_FISHE,
               CorrectedReads_MIFISHU,NormCorrectedReads_MIFISHU,CorrectedPercentReads_MIFISHU,NormPercReads_MIFISHU,NormPropReads_MIFISHU) %>% 
  pivot_longer(cols=c("CorrectedReads_12Steleo","NormCorrectedReads_12Steleo","CorrectedPercentReads_12Steleo","NormPercReads_12Steleo","NormPropReads_12Steleo",
                      "CorrectedReads_FISHE","NormCorrectedReads_FISHE","CorrectedPercentReads_FISHE","NormPercReads_FISHE","NormPropReads_FISHE",
                      "CorrectedReads_MIFISHU","NormCorrectedReads_MIFISHU","CorrectedPercentReads_MIFISHU","NormPercReads_MIFISHU","NormPropReads_MIFISHU"),
               names_to="DataType",
               values_to="Value") %>% 
  separate(col="DataType",
           into=c("DataType","Marker"),
           sep="_") %>% 
  mutate(DataType=factor(DataType),
         Marker=factor(Marker)) %>% 
  pivot_wider(names_from = "DataType",
              values_from = "Value") %>% 
  rename(CorrectedReadsPerLitre=NormCorrectedReads,
         PercentTotalReads=CorrectedPercentReads,
         PercentTotalReadsPerLitre=NormPercReads,
         ProportionTotalReadsPerLitre=NormPropReads) %>% 
  filter(CorrectedReads > 0) #remove any Sample/Marker combo where the Value is 0 (no reads)
         

##### AllMarkers Data plots and Outlier Removal#####
#picking data transformations
# ProportionTotalReadsPerLitre vs. NormQuantMean = P1
# CorrectedReadsPerLitre vs. NormQuantMean = P2
# ProportionTotalReadsPerLitre vs. log10(NormQuantMean) = P3
# CorrectedReadsPerLitre vs. log10(NormQuantMean) = P4
# log10(ProportionTotalReadsPerLitre) vs. log10(NormQuantMean) = P5
# log10(CorrectedReadsPerLitre) vs. log10(NormQuantMean) = P6


P1_base <- ggplot(AtlSalmon_data_long,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1 <- P1_base +
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
        panel.grid.major=element_blank());P1 #changes the axes, etc

P2_base <- ggplot(AtlSalmon_data_long,aes(x=NormQuantMean,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2 <- P2_base +
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
        panel.grid.major=element_blank());P2 #changes the axes, etc

P3_base <- ggplot(AtlSalmon_data_long,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3 <- P3_base +
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
        panel.grid.major=element_blank());P3 #changes the axes, etc

P4_base <- ggplot(AtlSalmon_data_long,aes(x=log10(NormQuantMean),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4 <- P4_base +
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
        panel.grid.major=element_blank());P4 #changes the axes, etc

P5_base <- ggplot(AtlSalmon_data_long,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5 <- P5_base +
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
        panel.grid.major=element_blank());P5 #changes the axes, etc

P6_base <- ggplot(AtlSalmon_data_long,aes(x=log10(NormQuantMean),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6 <- P6_base +
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
        panel.grid.major=element_blank());P6 #changes the axes, etc

#Save all 6 plots in 1
Multiplot <- grid.arrange(grobs=list(P1,P2,P3,P4,P5,P6),
                          cols=2,
                          top="Raw Data")
ggsave(Multiplot,
       file="DataVisualization_Transformations_AllMarkers_18Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_Prop <- AtlSalmon_data_long %>%  
  dplyr::select(!c(CorrectedReads,CorrectedReadsPerLitre,PercentTotalReads,PercentTotalReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_Prop$ProportionTotalReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_Prop$ProportionTotalReadsPerLitre),scientific=F)
# [1] "0.000001240196" "1.057051282051"
boxplot(AtlSalmon_data_Prop$ProportionTotalReadsPerLitre)
summary(AtlSalmon_data_Prop$ProportionTotalReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000012 0.1070659 0.2391554 0.2896911 0.4325539 1.0570513 
IQR <- IQR(AtlSalmon_data_Prop$ProportionTotalReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_Prop$ProportionTotalReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_NoOutliers_Temp <- AtlSalmon_data_Prop %>% 
  filter(ProportionTotalReadsPerLitre > Lower & ProportionTotalReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_NoOutliers_Temp$ProportionTotalReadsPerLitre) 
boxplot(AtlSalmon_data_Prop_NoOutliers_Temp$ProportionTotalReadsPerLitre)
#6 no outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_Prop_NoOutliers_Temp$NormQuantMean) #data is not normal, can't use z-score
range(AtlSalmon_data_Prop_NoOutliers_Temp$NormQuantMean)
# [1]   1.348685 248.628400
summary(AtlSalmon_data_Prop_NoOutliers_Temp$NormQuantMean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.441   7.991  12.952  17.857  21.131 248.628 
boxplot(AtlSalmon_data_Prop_NoOutliers_Temp$NormQuantMean)

IQR <- IQR(AtlSalmon_data_Prop_NoOutliers_Temp$NormQuantMean)
quartiles <- quantile(AtlSalmon_data_Prop_NoOutliers_Temp$NormQuantMean,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_NoOutliers <- AtlSalmon_data_Prop_NoOutliers_Temp %>% 
  filter(NormQuantMean > Lower & NormQuantMean < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_NoOutliers$NormQuantMean) 
boxplot(AtlSalmon_data_Prop_NoOutliers$NormQuantMean)
rm(AtlSalmon_data_Prop_NoOutliers_Temp)
#37 outliers removed


#redo plots with no outliers

P1_base_NoOut <- ggplot(AtlSalmon_data_Prop_NoOutliers,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_NoOut <- P1_base_NoOut +
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
        panel.grid.major=element_blank());P1 #changes the axes, etc

P3_base_NoOut <- ggplot(AtlSalmon_data_Prop_NoOutliers,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_NoOut <- P3_base_NoOut +
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
        panel.grid.major=element_blank());P3 #changes the axes, etc

P5_base_NoOut <- ggplot(AtlSalmon_data_Prop_NoOutliers,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_NoOut <- P5_base_NoOut +
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
        panel.grid.major=element_blank());P5 #changes the axes, etc



#Save all 6 plots in 1 (multiplot function from online)
Multiplot_NoOut <- grid.arrange(grobs=list(P1_NoOut,P1,P2,P3_NoOut,P3,P4,P5_NoOut,P5,P6),
                          cols=2,
                          top="No Outliers")
ggsave(Multiplot_NoOut,
       file="DataVisualization_Transformations_AllMarkers_NoOutliers_18Oct2023.pdf", 
       height=20, width=25,units = "in")





#####All Markers Error Structure Selection #####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

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
hist(AtlSalmon_data_Prop_NoOutliers$ProportionTotalReadsPerLitre) 
format(range(AtlSalmon_data_Prop_NoOutliers$ProportionTotalReadsPerLitre),scientific=F) #[1] "0.000001240196" "0.916276595745"
var(AtlSalmon_data_Prop_NoOutliers$ProportionTotalReadsPerLitre)# 0.04414822
mean(AtlSalmon_data_Prop_NoOutliers$ProportionTotalReadsPerLitre)# 0.2712005

descdist(AtlSalmon_data_Prop_NoOutliers$ProportionTotalReadsPerLitre, boot=500) #likely closest to normal or beta dist

#Normal error distribution
Full_LMM_normal_iden <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                           ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                           na.action = na.omit,
                           family =gaussian, #linear
                           REML=F)
Full_LMM_normal_log <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                           ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                           na.action = na.omit,
                           family =gaussian(link="log"), #logarithmic?
                           REML=F)
Full_LMM_normal_logit <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                               ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                               na.action = na.omit,
                               family =gaussian(link="logit"), #logarithmic?
                               REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the Full model with different error structures, then compare with diagnostic plots

#beta - fall between 0 and 1 (proportion!) 
Full_LMM_Beta <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                         ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                         na.action = na.omit,
                         family = beta_family()) #default logit

#Gamma
Full_LMM_Gamma <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                        ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                        family=Gamma)


#Model comparisons
#Plots - res vs fit and qqnorm
plot(simulateResiduals(Full_LMM_normal_iden))
  plot(x=fitted(Full_LMM_normal_iden),y=residuals(Full_LMM_normal_iden))
plot(simulateResiduals(Full_LMM_normal_log))
  plot(x=fitted(Full_LMM_normal_log),y=residuals(Full_LMM_normal_log))
plot(simulateResiduals(Full_LMM_normal_logit))
  plot(x=fitted(Full_LMM_normal_logit),y=residuals(Full_LMM_normal_logit))
plot(simulateResiduals(Full_LMM_Beta))
  plot(x=fitted(Full_LMM_Beta),y=residuals(Full_LMM_Beta))
plot(simulateResiduals(Full_LMM_Gamma)) #not great
  plot(x=fitted(Full_LMM_Gamma),y=residuals(Full_LMM_Gamma))



#histogram of residuals
hist(residuals(Full_LMM_normal_iden),main="Normal")
hist(residuals(Full_LMM_normal_log),main="Normal w/log")
hist(residuals(Full_LMM_normal_logit),main="Normal w/logit")
hist(residuals(Full_LMM_Beta),main="Beta")
hist(residuals(Full_LMM_Gamma),main="Gamma")

#test dispersion
testDispersion(simulateResiduals(Full_LMM_normal_iden,1000))
testDispersion(simulateResiduals(Full_LMM_normal_log,1000))
testDispersion(simulateResiduals(Full_LMM_normal_logit,1000))
testDispersion(simulateResiduals(Full_LMM_Beta,1000))
testDispersion(simulateResiduals(Full_LMM_Gamma,1000))

#AIC
AIC(Full_LMM_normal_iden)
AIC(Full_LMM_normal_log)
AIC(Full_LMM_normal_logit)
AIC(Full_LMM_Beta)
AIC(Full_LMM_Gamma)
#given the data, move forward with Beta into variable selection


#####All Markers Variable Selection#####

#full model
Full_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                    ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                    na.action = na.omit,
                    family = beta_family())
summary(Full_LMM)
AIC(Full_LMM) #-788.0444

#sort out random using AIC
Random1_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                       ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|RiverCode),
                       na.action = na.omit,
                       family = beta_family())
Random2_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                       ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|Run),
                       na.action = na.omit,
                       family = beta_family())
Random3_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, 
                       ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type,
                       na.action = na.omit,
                       family = beta_family())

#3 random models are: Full_LMM, Random1_LMM, Random2_LMM
summary(Full_LMM)    
# AIC      BIC   logLik deviance df.resid 
# -788.0   -730.1    407.0   -814.0      622  
summary(Random1_LMM) 
# AIC      BIC   logLik deviance df.resid 
# -790.0   -736.6    407.0   -814.0      623 
summary(Random2_LMM)    
# AIC      BIC   logLik deviance df.resid 
# -670.6   -617.2    347.3   -694.6      623   
summary(Random3_LMM)
# AIC      BIC   logLik deviance df.resid 
# -660.8   -611.8    341.4   -682.8      624   

#Random1_LMM has the best AIC score

#Select Fixed Effects, use AIC
Full_LMM_Fixed <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.fail,family = beta_family(),
                          ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|RiverCode))
AICc(Full_LMM_Fixed)

#use dredge to run all possible sub-models
Dredge_Full_LMM_Fixed <- dredge(Full_LMM_Fixed)
print(Dredge_Full_LMM_Fixed)
#Model24 is the best -> Log10 DNA, Marker, NQM, Year
#since the deltaAIC is so small between Model 24 and Model 23, go with the simpler model = Model 23 Marker, NQM, Year

##### AllMarkers Final model#####
FinalModel <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.fail,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + Year + (1|RiverCode))
summary(FinalModel) 
AICc(FinalModel)# -792.1818
Anova(FinalModel) 

#plain LM to compare
LM_finalmodel <- lm(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,
                    ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + Year + RiverCode)
summary(LM_finalmodel) 
AICc(LM_finalmodel) #-425.3581

#Final model with log normal
FinalModel_log <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = gaussian(link="logit"),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + Year + (1|RiverCode))
summary(FinalModel_log) 
AICc(FinalModel_log)# -421.3095

###plot the model
# Extract the prediction data frame
pred.mm <- ggpredict(FinalModel, terms = c("NormQuantMean"))  # this gives overall predictions for the model
pred.lm <- ggpredict(LM_finalmodel, terms = c("NormQuantMean"))
pred.log <- ggpredict(FinalModel_log, terms = c("NormQuantMean"))


# Plot the predictions 
FinalModel_plot <- ggplot(pred.mm) + 
  geom_line(data = pred.mm, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM\nProportionTotalReadsPerLitre ~ NormQuantMean + Marker + Year + (1|RiverCode)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalModel_plot

FinalLModel_plot <- ggplot(pred.lm) + 
  geom_line(data = pred.lm, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.lm,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total 12Steleo Reads\nPer Litre Filtered",
       title = "LM\nProportionTotalReadsPerLitre ~ NormQuantMean + Marker + Year + RiverCode") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLModel_plot

FinalLogModel_plot <- ggplot(pred.log) + 
  geom_line(data = pred.log, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.log,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total 12Steleo Reads\nPer Litre Filtered",
       title = "Log Normal MM\nProportionTotalReadsPerLitre ~ NormQuantMean + Marker + Year + (1|RiverCode)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLogModel_plot

#export plots
ggsave(FinalModel_plot, #plot you want to save
       file = "AtlanticSalmon_AllMarkers_GLMMBeta_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLModel_plot, #plot you want to save
       file = "AtlanticSalmon_AllMarkers_LM_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLogModel_plot, #plot you want to save
       file = "AtlanticSalmon_AllMarkers_LogNormalMM_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

########## 12S Data plots and Outlier Removal##########
#picking data transformations
# ProportionTotalReadsPerLitre vs. NormQuantMean = P1
# CorrectedReadsPerLitre vs. NormQuantMean = P2
# ProportionTotalReadsPerLitre vs. log10(NormQuantMean) = P3
# CorrectedReadsPerLitre vs. log10(NormQuantMean) = P4
# log10(ProportionTotalReadsPerLitre) vs. log10(NormQuantMean) = P5
# log10(CorrectedReadsPerLitre) vs. log10(NormQuantMean) = P6

AtlSalmon_data_long_12S <- AtlSalmon_data_long %>% 
  filter(AtlSalmon_data_long$Marker == "12Steleo")

P1_base_12S <- ggplot(AtlSalmon_data_long_12S,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_12S <- P1_base_12S +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P2_base_12S <- ggplot(AtlSalmon_data_long_12S,aes(x=NormQuantMean,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_12S <- P2_base_12S +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P3_base_12S <- ggplot(AtlSalmon_data_long_12S,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_12S <- P3_base_12S +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P4_base_12S <- ggplot(AtlSalmon_data_long_12S,aes(x=log10(NormQuantMean),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_12S <- P4_base_12S +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P5_base_12S <- ggplot(AtlSalmon_data_long_12S,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_12S <- P5_base_12S +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P6_base_12S <- ggplot(AtlSalmon_data_long_12S,aes(x=log10(NormQuantMean),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_12S <- P6_base_12S +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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
Multiplot_12S <- grid.arrange(grobs=list(P1_12S,P2_12S,P3_12S,P4_12S,P5_12S,P6_12S),
                          cols=2,
                          top="Raw Data 12S")
ggsave(Multiplot_12S,
       file="DataVisualization_Transformations_12S_18Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_Prop_12S <- AtlSalmon_data_long_12S %>%  
  dplyr::select(!c(CorrectedReads,CorrectedReadsPerLitre,PercentTotalReads,PercentTotalReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_Prop_12S$ProportionTotalReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_Prop_12S$ProportionTotalReadsPerLitre),scientific=F)
# [1] "0.000001365306" "1.057051282051"
boxplot(AtlSalmon_data_Prop_12S$ProportionTotalReadsPerLitre)
summary(AtlSalmon_data_Prop_12S$ProportionTotalReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000014 0.1105140 0.2692557 0.3023095 0.4307826 1.0570513 
IQR <- IQR(AtlSalmon_data_Prop_12S$ProportionTotalReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_Prop_12S$ProportionTotalReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_12S_NoOutliers_Temp <- AtlSalmon_data_Prop_12S %>% 
  filter(ProportionTotalReadsPerLitre > Lower & ProportionTotalReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_12S_NoOutliers_Temp$ProportionTotalReadsPerLitre) 
boxplot(AtlSalmon_data_Prop_12S_NoOutliers_Temp$ProportionTotalReadsPerLitre)
#2 no outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_Prop_12S_NoOutliers_Temp$NormQuantMean) #data is not normal, can't use z-score
range(AtlSalmon_data_Prop_12S_NoOutliers_Temp$NormQuantMean)
# [1]    1.441469 248.628400
summary(AtlSalmon_data_Prop_12S_NoOutliers_Temp$NormQuantMean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.441   7.851  12.609  18.275  20.950 248.628 
boxplot(AtlSalmon_data_Prop_12S_NoOutliers_Temp$NormQuantMean)

IQR <- IQR(AtlSalmon_data_Prop_12S_NoOutliers_Temp$NormQuantMean)
quartiles <- quantile(AtlSalmon_data_Prop_12S_NoOutliers_Temp$NormQuantMean,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_12S_NoOutliers <- AtlSalmon_data_Prop_12S_NoOutliers_Temp %>% 
  filter(NormQuantMean > Lower & NormQuantMean < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_12S_NoOutliers$NormQuantMean) 
boxplot(AtlSalmon_data_Prop_12S_NoOutliers$NormQuantMean)
rm(AtlSalmon_data_Prop_12S_NoOutliers_Temp)
#15 outliers removed


#redo plots with no outliers

P1_base_12S_NoOut <- ggplot(AtlSalmon_data_Prop_12S_NoOutliers,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_12S_NoOut <- P1_base_12S_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P3_base_12S_NoOut <- ggplot(AtlSalmon_data_Prop_12S_NoOutliers,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_12S_NoOut <- P3_base_12S_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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

P5_base_12S_NoOut <- ggplot(AtlSalmon_data_Prop_12S_NoOutliers,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_12S_NoOut <- P5_base_12S_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
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
Multiplot_12S_NoOut <- grid.arrange(grobs=list(P1_12S_NoOut,P1_12S,P2_12S,P3_12S_NoOut,P3_12S,P4_12S,P5_12S_NoOut,P5_12S,P6_12S),
                                cols=2,
                                top="No Outliers 12S")
ggsave(Multiplot_12S_NoOut,
       file="DataVisualization_Transformations_12S_NoOutliers_18Oct2023.pdf", 
       height=20, width=25,units = "in")


########## 12S GLM Variable Selection ##########
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation


#full model
Full_LMM_12S <- glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, 
                    ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                    na.action = na.omit,
                    family = beta_family())
summary(Full_LMM_12S)
AIC(Full_LMM_12S) #-284.4276

#sort out random using AIC
Random1_LMM_12S <- glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, 
                       ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                       na.action = na.omit,
                       family = beta_family())
Random2_LMM_12S <- glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, 
                       ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run),
                       na.action = na.omit,
                       family = beta_family())
Random3_LMM_12S <- glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, 
                       ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type,
                       na.action = na.omit,
                       family = beta_family())

#3 random models are: Full_LMM, Random1_LMM, Random2_LMM
summary(Full_LMM_12S)    
# AIC      BIC   logLik deviance df.resid 
# -284.4   -245.6    153.2   -306.4      242 
summary(Random1_LMM_12S) 
# AIC      BIC   logLik deviance df.resid 
# -286.4   -251.1    153.2   -306.4      243 
summary(Random2_LMM_12S)    
# AIC      BIC   logLik deviance df.resid 
# -241.9   -206.5    130.9   -261.9      243   
summary(Random3_LMM_12S)
# AIC      BIC   logLik deviance df.resid 
# -243.3   -211.5    130.7   -261.3      244 

#Random1_LMM has the best AIC score

#Select Fixed Effects, use AIC
Full_LMM_Fixed_12S <-  glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, 
                           ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                           na.action = na.fail,
                           family = beta_family())
AICc(Full_LMM_Fixed_12S) #-285.5186

#use dredge to run all possible sub-models
Dredge_Full_LMM_Fixed_12S <- dredge(Full_LMM_Fixed_12S)
print(Dredge_Full_LMM_Fixed_12S)
#Model 11 is the best -> NQM, Year



######### 12S Final model ##########
FinalModel_12S <- glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Year + (1|RiverCode))
summary(FinalModel_12S) 
AICc(FinalModel_12S)# -293.6206
Anova(FinalModel_12S)

#plain LM to compare
LM_finalmodel_12S <- lm(data = AtlSalmon_data_Prop_12S_NoOutliers, na.action = na.omit,
                    ProportionTotalReadsPerLitre ~ NormQuantMean + Year + RiverCode)
summary(LM_finalmodel_12S) 
AIC(LM_finalmodel_12S) #-178.6955

#Final model with log normal
FinalModel_12S_log <- glmmTMB(data = AtlSalmon_data_Prop_12S_NoOutliers, na.action = na.omit,family = gaussian(link="logit"),
                          ProportionTotalReadsPerLitre ~ NormQuantMean + Year + (1|RiverCode))
summary(FinalModel_12S_log) 
AIC(FinalModel_12S_log)# -156.7843

###plot the model
# Extract the prediction data frame
pred.mm_12S <- ggpredict(FinalModel_12S, terms = c("NormQuantMean"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_12S_plot <- ggplot(pred.mm_12S) + 
  geom_line(data = pred.mm, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_12S_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM 12S\nProportionTotalReadsPerLitre ~ NormQuantMean + Year + (1|RiverCode)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_12S_plot, #plot you want to save
       file = "AtlanticSalmon_12S_GLMMBeta_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


########## FISHE Data plots and Outlier Removal##########
#picking data transformations
# ProportionTotalReadsPerLitre vs. NormQuantMean = P1
# CorrectedReadsPerLitre vs. NormQuantMean = P2
# ProportionTotalReadsPerLitre vs. log10(NormQuantMean) = P3
# CorrectedReadsPerLitre vs. log10(NormQuantMean) = P4
# log10(ProportionTotalReadsPerLitre) vs. log10(NormQuantMean) = P5
# log10(CorrectedReadsPerLitre) vs. log10(NormQuantMean) = P6

AtlSalmon_data_long_FISHE <- AtlSalmon_data_long %>% 
  filter(AtlSalmon_data_long$Marker == "FISHE")

P1_base_FISHE <- ggplot(AtlSalmon_data_long_FISHE,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE <- P1_base_FISHE +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P2_base_FISHE <- ggplot(AtlSalmon_data_long_FISHE,aes(x=NormQuantMean,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_FISHE <- P2_base_FISHE +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P3_base_FISHE <- ggplot(AtlSalmon_data_long_FISHE,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE <- P3_base_FISHE +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P4_base_FISHE <- ggplot(AtlSalmon_data_long_FISHE,aes(x=log10(NormQuantMean),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_FISHE <- P4_base_FISHE +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P5_base_FISHE <- ggplot(AtlSalmon_data_long_FISHE,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE <- P5_base_FISHE +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P6_base_FISHE <- ggplot(AtlSalmon_data_long_FISHE,aes(x=log10(NormQuantMean),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_FISHE <- P6_base_FISHE +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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
Multiplot_FISHE <- grid.arrange(grobs=list(P1_FISHE,P2_FISHE,P3_FISHE,P4_FISHE,P5_FISHE,P6_FISHE),
                              cols=2,
                              top="Raw Data FISHE")
ggsave(Multiplot_FISHE,
       file="DataVisualization_Transformations_FISHE_18Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_Prop_FISHE <- AtlSalmon_data_long_FISHE %>%  
  dplyr::select(!c(CorrectedReads,CorrectedReadsPerLitre,PercentTotalReads,PercentTotalReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_Prop_FISHE$ProportionTotalReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_Prop_FISHE$ProportionTotalReadsPerLitre),scientific=F)
# [1] "0.005157143" "0.769381443"
boxplot(AtlSalmon_data_Prop_FISHE$ProportionTotalReadsPerLitre)
summary(AtlSalmon_data_Prop_FISHE$ProportionTotalReadsPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005157 0.102041 0.211919 0.268542 0.392200 0.769381 
IQR <- IQR(AtlSalmon_data_Prop_FISHE$ProportionTotalReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_Prop_FISHE$ProportionTotalReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_FISHE_NoOutliers_Temp <- AtlSalmon_data_Prop_FISHE %>% 
  filter(ProportionTotalReadsPerLitre > Lower & ProportionTotalReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$ProportionTotalReadsPerLitre) 
boxplot(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$ProportionTotalReadsPerLitre)
#no outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$NormQuantMean) #data is not normal, can't use z-score
range(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$NormQuantMean)
# [1]  1.906399 83.533156
summary(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$NormQuantMean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.906   8.757  13.932  17.217  21.428  83.533 
boxplot(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$NormQuantMean)

IQR <- IQR(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$NormQuantMean)
quartiles <- quantile(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp$NormQuantMean,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_FISHE_NoOutliers <- AtlSalmon_data_Prop_FISHE_NoOutliers_Temp %>% 
  filter(NormQuantMean > Lower & NormQuantMean < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_FISHE_NoOutliers$NormQuantMean) 
boxplot(AtlSalmon_data_Prop_FISHE_NoOutliers$NormQuantMean)
rm(AtlSalmon_data_Prop_FISHE_NoOutliers_Temp)
#10 outliers removed


#redo plots with no outliers

P1_base_FISHE_NoOut <- ggplot(AtlSalmon_data_Prop_FISHE_NoOutliers,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoOut <- P1_base_FISHE_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P3_base_FISHE_NoOut <- ggplot(AtlSalmon_data_Prop_FISHE_NoOutliers,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoOut <- P3_base_FISHE_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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

P5_base_FISHE_NoOut <- ggplot(AtlSalmon_data_Prop_FISHE_NoOutliers,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoOut <- P5_base_FISHE_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
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
Multiplot_FISHE_NoOut <- grid.arrange(grobs=list(P1_FISHE_NoOut,P1_FISHE,P2_FISHE,P3_FISHE_NoOut,P3_FISHE,P4_FISHE,P5_FISHE_NoOut,P5_FISHE,P6_FISHE),
                                    cols=2,
                                    top="No Outliers FISHE")
ggsave(Multiplot_FISHE_NoOut,
       file="DataVisualization_Transformations_FISHE_NoOutliers_18Oct2023.pdf", 
       height=20, width=25,units = "in")


########## FISHE GLM Variable Selection ##########
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#full model
Full_LMM_FISHE <- glmmTMB(data = AtlSalmon_data_Prop_FISHE_NoOutliers, 
                        ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                        na.action = na.omit,
                        family = beta_family())
summary(Full_LMM_FISHE)
AICc(Full_LMM_FISHE) #-190.5722

#sort out random using AIC
Random1_LMM_FISHE <- glmmTMB(data = AtlSalmon_data_Prop_FISHE_NoOutliers, 
                           ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                           na.action = na.omit,
                           family = beta_family())
Random2_LMM_FISHE <- glmmTMB(data = AtlSalmon_data_Prop_FISHE_NoOutliers, 
                           ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run),
                           na.action = na.omit,
                           family = beta_family())
Random3_LMM_FISHE <- glmmTMB(data = AtlSalmon_data_Prop_FISHE_NoOutliers, 
                           ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type,
                           na.action = na.omit,
                           family = beta_family())

#3 random models are: Full_LMM, Random1_LMM, Random2_LMM
summary(Full_LMM_FISHE)    
# AIC      BIC   logLik deviance df.resid
# -192.3   -158.0    107.1   -214.3      156
summary(Random1_LMM_FISHE) 
# AIC      BIC   logLik deviance df.resid
# -194.3   -163.1    107.1   -214.3      157 
summary(Random2_LMM_FISHE)    
# AIC      BIC   logLik deviance df.resid 
# -187.1   -155.9    103.5   -207.1      157   
summary(Random3_LMM_FISHE)
# AIC      BIC   logLik deviance df.resid 
# -188.9   -160.9    103.5   -206.9      158 

#Random1_LMM has the best AIC score

#Select Fixed Effects, use AIC
Full_LMM_Fixed_FISHE <-  glmmTMB(data = AtlSalmon_data_Prop_FISHE_NoOutliers, 
                               ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                               na.action = na.fail,
                               family = beta_family())
AICc(Full_LMM_Fixed_FISHE) #-192.8651

#use dredge to run all possible sub-models
Dredge_Full_LMM_Fixed_FISHE <- dredge(Full_LMM_Fixed_FISHE)
print(Dredge_Full_LMM_Fixed_FISHE)
#Model 12 is the best -> NQM, Year, log10DNAConc



######### FISHE Final model ##########
FinalModel_FISHE <- glmmTMB(data = AtlSalmon_data_Prop_FISHE_NoOutliers, na.action = na.omit,family = beta_family(),
                          ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + (1|RiverCode))
summary(FinalModel_FISHE) 
AICc(FinalModel_FISHE)# -193.1442
Anova(FinalModel_FISHE)

###plot the model
# Extract the prediction data frame
pred.mm_FISHE <- ggpredict(FinalModel_FISHE, terms = c("NormQuantMean"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_FISHE_plot <- ggplot(pred.mm_FISHE) + 
  geom_line(data = pred.mm, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_FISHE_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM FISHE\nProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + (1|RiverCode)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_FISHE_plot, #plot you want to save
       file = "AtlanticSalmon_FISHE_GLMMBeta_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

########## MIFISHU Data plots and Outlier Removal##########
#picking data transformations
# ProportionTotalReadsPerLitre vs. NormQuantMean = P1
# CorrectedReadsPerLitre vs. NormQuantMean = P2
# ProportionTotalReadsPerLitre vs. log10(NormQuantMean) = P3
# CorrectedReadsPerLitre vs. log10(NormQuantMean) = P4
# log10(ProportionTotalReadsPerLitre) vs. log10(NormQuantMean) = P5
# log10(CorrectedReadsPerLitre) vs. log10(NormQuantMean) = P6

AtlSalmon_data_long_MIFISHU <- AtlSalmon_data_long %>% 
  filter(AtlSalmon_data_long$Marker == "MIFISHU")

P1_base_MIFISHU <- ggplot(AtlSalmon_data_long_MIFISHU,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU <- P1_base_MIFISHU +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P2_base_MIFISHU <- ggplot(AtlSalmon_data_long_MIFISHU,aes(x=NormQuantMean,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_MIFISHU <- P2_base_MIFISHU +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P3_base_MIFISHU <- ggplot(AtlSalmon_data_long_MIFISHU,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU <- P3_base_MIFISHU +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P4_base_MIFISHU <- ggplot(AtlSalmon_data_long_MIFISHU,aes(x=log10(NormQuantMean),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_MIFISHU <- P4_base_MIFISHU +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P5_base_MIFISHU <- ggplot(AtlSalmon_data_long_MIFISHU,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU <- P5_base_MIFISHU +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P6_base_MIFISHU <- ggplot(AtlSalmon_data_long_MIFISHU,aes(x=log10(NormQuantMean),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_MIFISHU <- P6_base_MIFISHU +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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
Multiplot_MIFISHU <- grid.arrange(grobs=list(P1_MIFISHU,P2_MIFISHU,P3_MIFISHU,P4_MIFISHU,P5_MIFISHU,P6_MIFISHU),
                                cols=2,
                                top="Raw Data MIFISHU")
ggsave(Multiplot_MIFISHU,
       file="DataVisualization_Transformations_MIFISHU_18Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_Prop_MIFISHU <- AtlSalmon_data_long_MIFISHU %>%  
  dplyr::select(!c(CorrectedReads,CorrectedReadsPerLitre,PercentTotalReads,PercentTotalReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_Prop_MIFISHU$ProportionTotalReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_Prop_MIFISHU$ProportionTotalReadsPerLitre),scientific=F)
# [1] "0.000001240196" "0.967790697674"
boxplot(AtlSalmon_data_Prop_MIFISHU$ProportionTotalReadsPerLitre)
summary(AtlSalmon_data_Prop_MIFISHU$ProportionTotalReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000012 0.1038000 0.2344898 0.2911476 0.4472046 0.9677907
IQR <- IQR(AtlSalmon_data_Prop_MIFISHU$ProportionTotalReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_Prop_MIFISHU$ProportionTotalReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp <- AtlSalmon_data_Prop_MIFISHU %>% 
  filter(ProportionTotalReadsPerLitre > Lower & ProportionTotalReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$ProportionTotalReadsPerLitre) 
boxplot(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$ProportionTotalReadsPerLitre)
#1 outlier removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$NormQuantMean) #data is not normal, can't use z-score
range(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$NormQuantMean)
# [1]  1.441469 248.628400
summary(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$NormQuantMean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.906   8.757  13.932  17.217  21.428  83.533 
boxplot(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$NormQuantMean)

IQR <- IQR(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$NormQuantMean)
quartiles <- quantile(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp$NormQuantMean,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_Prop_MIFISHU_NoOutliers <- AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp %>% 
  filter(NormQuantMean > Lower & NormQuantMean < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_Prop_MIFISHU_NoOutliers$NormQuantMean) 
boxplot(AtlSalmon_data_Prop_MIFISHU_NoOutliers$NormQuantMean)
rm(AtlSalmon_data_Prop_MIFISHU_NoOutliers_Temp)
#12 outliers removed


#redo plots with no outliers

P1_base_MIFISHU_NoOut <- ggplot(AtlSalmon_data_Prop_MIFISHU_NoOutliers,aes(x=NormQuantMean,y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoOut <- P1_base_MIFISHU_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P3_base_MIFISHU_NoOut <- ggplot(AtlSalmon_data_Prop_MIFISHU_NoOutliers,aes(x=log10(NormQuantMean),y=ProportionTotalReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoOut <- P3_base_MIFISHU_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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

P5_base_MIFISHU_NoOut <- ggplot(AtlSalmon_data_Prop_MIFISHU_NoOutliers,aes(x=log10(NormQuantMean),y=log10(ProportionTotalReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoOut <- P5_base_MIFISHU_NoOut +
  geom_point(alpha=0.7) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
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
Multiplot_MIFISHU_NoOut <- grid.arrange(grobs=list(P1_MIFISHU_NoOut,P1_MIFISHU,P2_MIFISHU,P3_MIFISHU_NoOut,P3_MIFISHU,P4_MIFISHU,P5_MIFISHU_NoOut,P5_MIFISHU,P6_MIFISHU),
                                      cols=2,
                                      top="No Outliers MIFISHU")
ggsave(Multiplot_MIFISHU_NoOut,
       file="DataVisualization_Transformations_MIFISHU_NoOutliers_18Oct2023.pdf", 
       height=20, width=20,units = "in")


########## MIFISHU GLM Variable Selection ##########
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation


#full model
Full_LMM_MIFISHU <- glmmTMB(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, 
                          ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                          na.action = na.omit,
                          family = beta_family())
summary(Full_LMM_MIFISHU)
AICc(Full_LMM_MIFISHU) #-288.0136

#sort out random using AIC
Random1_LMM_MIFISHU <- glmmTMB(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, 
                             ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                             na.action = na.omit,
                             family = beta_family())
Random2_LMM_MIFISHU <- glmmTMB(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, 
                             ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run),
                             na.action = na.omit,
                             family = beta_family())
Random3_LMM_MIFISHU <- glmmTMB(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, 
                             ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type,
                             na.action = na.omit,
                             family = beta_family())

#3 random models are: Full_LMM, Random1_LMM, Random2_LMM
summary(Full_LMM_MIFISHU)    
# AIC      BIC   logLik deviance df.resid 
# -289.3   -252.1    155.6   -311.3      207 
summary(Random1_LMM_MIFISHU) 
# AIC      BIC   logLik deviance df.resid 
# -291.3   -257.5    155.6   -311.3      208 
summary(Random2_LMM_MIFISHU)    
# AIC      BIC   logLik deviance df.resid 
# -227.0   -193.1    123.5   -247.0      208    
summary(Random3_LMM_MIFISHU)
# AIC      BIC   logLik deviance df.resid 
# -228.7   -198.2    123.3   -246.7      209 

#Random1_LMM has the best AIC score

#Select Fixed Effects, use AIC
Full_LMM_Fixed_MIFISHU <-  glmmTMB(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, 
                                 ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                                 na.action = na.fail,
                                 family = beta_family())
AICc(Full_LMM_Fixed_MIFISHU) #-290.2324

#use dredge to run all possible sub-models
Dredge_Full_LMM_Fixed_MIFISHU <- dredge(Full_LMM_Fixed_MIFISHU)
print(Dredge_Full_LMM_Fixed_MIFISHU)
#Model 16 and 12 are the best, with minimal change between them -> NQM, Year, log10DNAConc, Type (only in Model 16)
#go with simpler model



######### MIFISHU Final model ##########
FinalModel_MIFISHU <- glmmTMB(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, na.action = na.omit,family = beta_family(),
                            ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + (1|RiverCode))
summary(FinalModel_MIFISHU) 
AICc(FinalModel_MIFISHU)# -289.9664
Anova(FinalModel_MIFISHU)

###plot the model
# Extract the prediction data frame
pred.mm_MIFISHU <- ggpredict(FinalModel_MIFISHU, terms = c("NormQuantMean"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_MIFISHU_plot <- ggplot(pred.mm_MIFISHU) + 
  geom_line(data = pred.mm, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_MIFISHU_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM MIFISHU\nProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + (1|RiverCode)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_MIFISHU_plot, #plot you want to save
       file = "AtlanticSalmon_MIFISHU_GLMMBeta_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


######Multiplot All FinalModel Plots#####
#Save all 6 plots in 1
Final_Multiplot <- grid.arrange(grobs=list(FinalModel_plot,
                                           FinalModel_12S_plot,
                                           FinalModel_FISHE_plot,
                                           FinalModel_MIFISHU_plot),
                          cols=2,
                          top="Beta GLMM")
ggsave(Final_Multiplot,
       file="DataVisualization_Transformations_GLMMMultiplot_18Oct2023.pdf", 
       height=20, width=25,units = "in")


#######Coefficient Estimates all models#######
#FinalModel
FinalModelSum <- summary(FinalModel)

lme4::formatVC(FinalModelSum$varcor$cond) #"random effects variances"
coef(FinalModelSum)$cond #"conditional fixed effects")

tab_model(FinalModel, show.ci = T)






FinalModel_12S


FinalModel_FISHE



FinalModel_MIFISHU



######save workspace######
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
