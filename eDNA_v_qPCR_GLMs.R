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
         

#####Data plots and Outlier Removal#####
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
       file="DataVisualization_Transformations_18Oct2023.pdf", 
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
  scale_fill_brewer(palette="Dark2")+
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
       file="DataVisualization_Transformations_NoOutliers_18Oct2023.pdf", 
       height=20, width=20,units = "in")





#####GLMs#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

AtlSalmon_data_Prop_NoOutliers


#following some stuff from Zuur 2009
#make dot charts of the continuous variables
op <- par(mfrow=c(3,2),mar=c(3,3,3,1))
dotchart(AtlSalmon_data_Prop_NoOutliers$ProportionTotalReadsPerLitre, main="NormPropReads",group=AtlSalmon_data_Prop_NoOutliers$RiverCode)
plot(0,0,type="n", axes=F)
dotchart(AtlSalmon_data_Prop_NoOutliers$NormQuantMean, main="MeanQuant",group=AtlSalmon_data_Prop_NoOutliers$RiverCode)
dotchart(log10(AtlSalmon_data_Prop_NoOutliers$NormQuantMean), main="Log10MeanQuant",group=AtlSalmon_data_Prop_NoOutliers$RiverCode)
dotchart(AtlSalmon_data_Prop_NoOutliers$DNAConcScale, main="Scaled DNAConc",group=AtlSalmon_data_Prop_NoOutliers$RiverCode)
dotchart(AtlSalmon_data_Prop_NoOutliers$log10DNAConc, main="Log10DNAConc",group=AtlSalmon_data_Prop_NoOutliers$RiverCode)
par(op)

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
Full_LMM_Fixed <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                          ProportionTotalReadsPerLitre ~ NormQuantMean + Marker + log10DNAConc + Year + Type + (1|RiverCode))


Fixed1_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Year + (1|RiverCode))
Fixed2_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + Type + (1|RiverCode))
Fixed3_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + log10DNAConc + (1|RiverCode))
Fixed4_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Year + Type + (1|RiverCode))
Fixed5_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Year + (1|RiverCode))
Fixed6_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Type + (1|RiverCode))
Fixed7_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + (1|RiverCode))
Fixed8_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ log10DNAConc + Year + Type + (1|RiverCode))
Fixed9_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ log10DNAConc + Year + (1|RiverCode))
Fixed10_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                       ProportionTotalReadsPerLitre ~ log10DNAConc + Type + (1|RiverCode))
Fixed11_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                       ProportionTotalReadsPerLitre ~ log10DNAConc + (1|RiverCode))
Fixed12_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                       ProportionTotalReadsPerLitre ~ Year + Type + (1|RiverCode))
Fixed13_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                       ProportionTotalReadsPerLitre ~ Year + (1|RiverCode))
Fixed14_LMM <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                       ProportionTotalReadsPerLitre ~ Type + (1|RiverCode))

#14 fixed models
summary(Full_LMM_Fixed) #AIC   -286.4   
summary(Fixed1_LMM) #AIC    -292.1    
summary(Fixed2_LMM) #AIC     -277.1   
summary(Fixed3_LMM) #AIC    -282.5   
summary(Fixed4_LMM) #AIC    -288.4   
summary(Fixed5_LMM) #AIC    -294.0   ***
summary(Fixed6_LMM) #AIC    -278.6   
summary(Fixed7_LMM) #AIC    -284.0   
summary(Fixed8_LMM) #AIC    -258.2   
summary(Fixed9_LMM) #AIC    -263.9    
summary(Fixed10_LMM) #AIC   -252.4   
summary(Fixed11_LMM) #AIC   -258.0
summary(Fixed12_LMM) #AIC   -257.7
summary(Fixed13_LMM) #AIC   -263.5
summary(Fixed14_LMM) #AIC   -248.8

#Note = LL will be better for larger models, AIC corrects for this
#best model by AIC selection is Fixed5_LMM, ProportionTotalReadsPerLitre ~ log10NormQuantMean + Year + (1|RiverCode)

#####Final model#####
FinalModel <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = beta_family(),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Year + (1|RiverCode))
summary(FinalModel) 
AIC(FinalModel)# -293.9621

#plain LM to compare
LM_finalmodel <- lm(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,
                    ProportionTotalReadsPerLitre ~ NormQuantMean + Year + RiverCode)
summary(LM_finalmodel) 
AIC(LM_finalmodel) #-178.6955

#Final model with log normal
FinalModel_log <- glmmTMB(data = AtlSalmon_data_Prop_NoOutliers, na.action = na.omit,family = gaussian(link="logit"),
                      ProportionTotalReadsPerLitre ~ NormQuantMean + Year + (1|RiverCode))
summary(FinalModel_log) 
AIC(FinalModel_log)# -156.7843

###plot the model
# Extract the prediction data frame
pred.mm <- ggpredict(FinalModel, terms = c("NormQuantMean"))  # this gives overall predictions for the model
pred.lm <- ggpredict(LM_finalmodel, terms = c("NormQuantMean"))
pred.log <- ggpredict(FinalModel_log, terms = c("NormQuantMean"))

# Plot the predictions 
FinalModel_plot <- ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_NoOutliers,                      # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre)) + 
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total 12Steleo Reads\nPer Litre Filtered",
       title = "Beta GLLM") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalModel_plot

FinalLModel_plot <- ggplot(pred.lm) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_NoOutliers,                      # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre)) + 
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total 12Steleo Reads\nPer Litre Filtered",
       title = "LM") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLModel_plot

FinalLogModel_plot <- ggplot(pred.log) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_Prop_NoOutliers,                      # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = ProportionTotalReadsPerLitre)) + 
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total 12Steleo Reads\nPer Litre Filtered",
       title = "Log Normal MM") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLogModel_plot

#export plots
ggsave(FinalModel_plot, #plot you want to save
       file = "AtlanticSalmon_12S_GLMMBeta_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLModel_plot, #plot you want to save
       file = "AtlanticSalmon_12S_LM_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

ggsave(FinalLogModel_plot, #plot you want to save
       file = "AtlanticSalmon_12S_LogNormalMM_18Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
