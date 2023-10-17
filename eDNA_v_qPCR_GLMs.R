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

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load data (see eDNA_v_qPCR_DataProcessing.R)####
AtlSalmon_data_raw <- read.csv("qPCR_eDNA_SalmoSalar_CleanedData_11Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

PinkSalmon_data <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_11Oct2023.csv",stringsAsFactors = T) %>% 
  dplyr::select(!(X))

ArcticCharr_data <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_11Oct2023.csv",stringsAsFactors = T) %>% 
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
         DNAConcScale=scale(DNAConc_pg_uL, center=T, scale=T),
         log10DNAConc=log10(DNAConc_pg_uL),
         log10NormQuantMean=log10(NormQuantMean)) 


AtlSalmon_data_12Steleo <- AtlSalmon_data %>% 
  dplyr::select(SampleID,
                Type,
                CorrectedReads_12Steleo,
                NormCorrectedReads_12Steleo,
                CorrectedPercentReads_12Steleo,
                NormPercReads_12Steleo,
                VolFiltered,
                DNAConc_pg_uL,
                DNAConcScale,
                log10DNAConc,
                Run,
                QuantMean,
                NormQuantMean,
                log10NormQuantMean,
                RiverCode,
                Year) %>% 
  filter(CorrectedReads_12Steleo > 0) %>%  #remove any samples that had no 12S reads
  filter(NormPercReads_12Steleo <= 100) %>%  #remove samples with >100% reads
  mutate(NormPropReads_12Steleo=NormPercReads_12Steleo/100) %>% 
  relocate(NormPropReads_12Steleo, .after = NormPercReads_12Steleo)

###remove outliers
#Response variable - proportion of total reads
hist(AtlSalmon_data_12Steleo$NormPropReads_12Steleo) #data is not normal, can't use z-score
range(AtlSalmon_data_12Steleo$NormPropReads_12Steleo)
#[1] 0.0000000 0.8662353
summary(AtlSalmon_data_12Steleo$NormPropReads_12Steleo)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.06701 0.21890 0.26078 0.41000 0.86624 
IQR <- IQR(AtlSalmon_data_12Steleo$NormPropReads_12Steleo)
quartiles <- quantile(AtlSalmon_data_12Steleo$NormPropReads_12Steleo,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_12Steleo_NoOutliers_Prep <- AtlSalmon_data_12Steleo %>% 
  filter(NormPropReads_12Steleo > Lower & NormPropReads_12Steleo < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormPropReads_12Steleo) 
boxplot(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormPropReads_12Steleo)
#no outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormQuantMean) #data is not normal, can't use z-score
range(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormQuantMean)
# [1]   1.348685 248.628400
summary(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormQuantMean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.349   6.848  11.626  17.633  19.523 248.628
boxplot(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormQuantMean)

IQR <- IQR(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormQuantMean)
quartiles <- quantile(AtlSalmon_data_12Steleo_NoOutliers_Prep$NormQuantMean,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_12Steleo_NoOutliers <- AtlSalmon_data_12Steleo_NoOutliers_Prep %>% 
  filter(NormQuantMean > Lower & NormQuantMean < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_12Steleo_NoOutliers$NormQuantMean) 
boxplot(AtlSalmon_data_12Steleo_NoOutliers$NormQuantMean)
rm(AtlSalmon_data_12Steleo_NoOutliers_Prep)
#15 outliers removed

#####Data plots#####

#picking data transformations
# NormPropReads_12Steleo vs. NormQuantMean = P1
# NormCorrectedReads_12Steleo vs. NormQuantMean = P2
# log10(NormPropReads_12Steleo) vs. log10(NormQuantMean) = P3
# log10(NormCorrectedReads_12Steleo) vs. log10(NormQuantMean) = P4
# log10(NormPropReads_12Steleo) vs. log10(NormQuantMean) = P5
# log10(NormCorrectedReads_12Steleo) vs. log10(NormQuantMean) = P6

P1_base <- ggplot(AtlSalmon_data_12Steleo_NoOutliers,aes(x=NormQuantMean,y=NormPropReads_12Steleo))
P1 <- P1_base +
  geom_point() +
  # scale_colour_manual(values = paletteer_d("ggthemes::Tableau_10"))+
  labs(x="Mean Copies\nPer Litre Filtered",
       y="Prop Total 12Steleo Reads\nPer Litre Filtered")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P1 #changes the axes, etc

P2_base <- ggplot(AtlSalmon_data_12Steleo_NoOutliers,aes(x=NormQuantMean,y=NormCorrectedReads_12Steleo))
P2 <- P2_base +
  geom_point() +
  labs(x="Mean Copies\nPer Litre Filtered",
       y="12Steleo Reads\nPer Litre Filtered")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P2 #changes the axes, etc

P3_base <- ggplot(AtlSalmon_data_12Steleo_NoOutliers,aes(x=log10(NormQuantMean),y=NormPropReads_12Steleo))
P3 <- P3_base +
  geom_point() +
  # scale_colour_manual(values = paletteer_d("ggthemes::Tableau_10"))+
    labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "Prop Total 12Steleo Reads\nPer Litre Filtered") +
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P3 #changes the axes, etc

P4_base <- ggplot(AtlSalmon_data_12Steleo_NoOutliers,aes(x=log10(NormQuantMean),y=NormCorrectedReads_12Steleo))
P4 <- P4_base +
  geom_point() +
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "12Steleo Reads\nPer Litre Filtered") +
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P4 #changes the axes, etc

P5_base <- ggplot(AtlSalmon_data_12Steleo_NoOutliers,aes(x=log10(NormQuantMean),y=log10(NormPropReads_12Steleo)))
P5 <- P5_base +
  geom_point() +
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(Prop Total 12Steleo Reads)", paste("Per Litre Filtered")))) +
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P5 #changes the axes, etc

P6_base <- ggplot(AtlSalmon_data_12Steleo_NoOutliers,aes(x=log10(NormQuantMean),y=log10(NormCorrectedReads_12Steleo)))
P6 <- P6_base +
  geom_point() +
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(12Steleo Reads)", paste("Per Litre Filtered")))) +
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P6 #changes the axes, etc

#Save all 6 plots in 1 (multiplot function from online)
dev.off()
pdf("DataVisualization_Transformations_17Oct2023.pdf", height=20, width=15)
multiplot(P1,P3,P5,P2,P4,P6,
          cols=2)
dev.off()

#####GLMs#####
#following https://ourcodingclub.github.io/tutorials/mixed-models/
###to start from Zuur 2009
# 1. fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# 2. sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# 3. sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# 4. once you arrive at the final model present it using REML estimation

#following some stuff from Zuur 2009
#make dot charts of the continuous variables
op <- par(mfrow=c(3,2),mar=c(3,3,3,1))
dotchart(AtlSalmon_data_12Steleo_NoOutliers$NormPropReads_12Steleo, main="NormPropReads",group=AtlSalmon_data_12Steleo_NoOutliers$RiverCode)
plot(0,0,type="n", axes=F)
dotchart(AtlSalmon_data_12Steleo_NoOutliers$NormQuantMean, main="MeanQuant",group=AtlSalmon_data_12Steleo_NoOutliers$RiverCode)
dotchart(AtlSalmon_data_12Steleo_NoOutliers$log10NormQuantMean, main="Log10MeanQuant",group=AtlSalmon_data_12Steleo_NoOutliers$RiverCode)
dotchart(AtlSalmon_data_12Steleo_NoOutliers$DNAConcScale, main="Scaled DNAConc",group=AtlSalmon_data_12Steleo_NoOutliers$RiverCode)
dotchart(AtlSalmon_data_12Steleo_NoOutliers$log10DNAConc, main="Log10DNAConc",group=AtlSalmon_data_12Steleo_NoOutliers$RiverCode)
par(op)

#Vars
# River = Random
# Type = Fixed - only 4 types, and should not be nested within river
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed

#picking error structure
#hist of response variable
hist(AtlSalmon_data_12Steleo_NoOutliers$NormPropReads_12Steleo) 
range(AtlSalmon_data_12Steleo_NoOutliers$NormPropReads_12Steleo)
var(AtlSalmon_data_12Steleo_NoOutliers$NormPropReads_12Steleo) #0.04398729
mean(AtlSalmon_data_12Steleo_NoOutliers$NormPropReads_12Steleo) #0.2827938

descdist(AtlSalmon_data_12Steleo_NoOutliers$NormPropReads_12Steleo, boot=500) #likely closest to normal or beta dist

#Normal error distribution
Full_LMM_normal_iden <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                           NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                           na.action = na.omit,
                           family =gaussian, #linear
                           REML=F)
Full_LMM_normal_log <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                           NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                           na.action = na.omit,
                           family =gaussian(link="log"), #logarithmic?
                           REML=F)
Full_LMM_normal_logit <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                               NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                               na.action = na.omit,
                               family =gaussian(link="logit"), #logarithmic?
                               REML=F)
#we have a bounded outcome in that the response variable is a proportion
#redo the Full model with different error structures, then compare with diagnostic plots

#beta - fall between 0 and 1 (proportion!) 
Full_LMM_Beta <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                         NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                         na.action = na.omit,
                         family = beta_family()) #default logit

#Gamma
Full_LMM_Gamma <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                        NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
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

AIC(Full_LMM_normal_log)
AIC(Full_LMM_Beta)
#given the data, move forward with Beta into variable selection

#full model
Full_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                    NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run) + (1|RiverCode),
                    na.action = na.omit,
                    family =beta_family())
summary(Full_LMM)
AIC(Full_LMM) #-284.4276

#sort out random using AIC
Random1_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                       NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode),
                       na.action = na.omit,
                       family = beta_family())
Random2_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                       NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|Run),
                       na.action = na.omit,
                       family = beta_family())
Random3_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, 
                       NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type,
                       na.action = na.omit,
                       family = beta_family())

#3 random models are: Full_LMM, Random1_LMM, Random2_LMM
summary(Full_LMM)    
# AIC      BIC   logLik deviance df.resid 
# -284.4   -245.6    153.2   -306.4      242  
summary(Random1_LMM) 
# AIC      BIC   logLik deviance df.resid 
# -286.4   -251.1    153.2   -306.4      243 
summary(Random2_LMM)    
# AIC      BIC   logLik deviance df.resid 
# -241.9   -206.5    130.9   -261.9      243  
summary(Random3_LMM)
# AIC      BIC   logLik deviance df.resid 
# -243.3   -211.5    130.7   -261.3      244  

#Random1_LMM has the best AIC score

#Select Fixed Effects, use AIC
Full_LMM_Fixed <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                          NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + Type + (1|RiverCode))
Fixed1_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Year + (1|RiverCode))
Fixed2_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + Type + (1|RiverCode))
Fixed3_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + log10DNAConc + (1|RiverCode))
Fixed4_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + Year + Type + (1|RiverCode))
Fixed5_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + Year + (1|RiverCode))
Fixed6_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + Type + (1|RiverCode))
Fixed7_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + (1|RiverCode))
Fixed8_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ log10DNAConc + Year + Type + (1|RiverCode))
Fixed9_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ log10DNAConc + Year + (1|RiverCode))
Fixed10_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                       NormPropReads_12Steleo ~ log10DNAConc + Type + (1|RiverCode))
Fixed11_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                       NormPropReads_12Steleo ~ log10DNAConc + (1|RiverCode))
Fixed12_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                       NormPropReads_12Steleo ~ Year + Type + (1|RiverCode))
Fixed13_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                       NormPropReads_12Steleo ~ Year + (1|RiverCode))
Fixed14_LMM <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                       NormPropReads_12Steleo ~ Type + (1|RiverCode))

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
#best model by AIC selection is Fixed5_LMM, NormPropReads_12Steleo ~ log10NormQuantMean + Year + (1|RiverCode)

###Final model
FinalModel <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = beta_family(),
                      NormPropReads_12Steleo ~ NormQuantMean + Year + (1|RiverCode))
summary(FinalModel) 
AIC(FinalModel)# -293.9621

#plain LM to compare
LM_finalmodel <- lm(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,
                    NormPropReads_12Steleo ~ NormQuantMean + Year + RiverCode)
summary(LM_finalmodel) 
AIC(LM_finalmodel) #-178.6955

#Final model with log normal
FinalModel_log <- glmmTMB(data = AtlSalmon_data_12Steleo_NoOutliers, na.action = na.omit,family = gaussian(link="logit"),
                      NormPropReads_12Steleo ~ NormQuantMean + Year + (1|RiverCode))
summary(FinalModel_log) 
AIC(FinalModel_log)# -293.9621

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
  geom_point(data = AtlSalmon_data_12Steleo_NoOutliers,                      # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = NormPropReads_12Steleo)) + 
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
        panel.grid.major=element_blank());FinalModel_plot

FinalLModel_plot <- ggplot(pred.lm) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_12Steleo_NoOutliers,                      # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = NormPropReads_12Steleo)) + 
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
  geom_point(data = AtlSalmon_data_12Steleo_NoOutliers,                      # adding the raw data (scaled values)
             aes(x = NormQuantMean, y = NormPropReads_12Steleo)) + 
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total 12Steleo Reads\nPer Litre Filtered",
       title = "Log Normal") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());FinalLogModel_plot




#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
