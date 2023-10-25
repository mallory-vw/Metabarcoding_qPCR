#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GAMs_workspace.RData")

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
library("mgcv")

#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load data (see eDNA_v_qPCR_DataProcessing.R and eDNA_v_qPCR_GLMs.R)####
###Zero read counts removed###
AllSpecies_AllMarkers_NoZero_NoOutliers <- read.csv("AllSpecies_AllMarkers_NoZero_NoOutliers_25Oct2023.csv", stringsAsFactors = T) %>% 
  dplyr::select(!(X))

AtlSalmon_AllMarkers_NoZero_NoOutliers <- read.csv("AtlSalmon_AllMarkers_NoZero_NoOutliers_25Oct2023.csv", stringsAsFactors = T) %>% 
  dplyr::select(!(X))

ArcCharr_AllMarkers_NoZero_NoOutliers <- read.csv("ArcticCharr_AllMarkers_NoZero_NoOutliers_25Oct2023.csv", stringsAsFactors = T) %>% 
  dplyr::select(!(X))

PinkSalmon_AllMarkers_NoZero_NoOutliers <- read.csv("PinkSalmon_AllMarkers_NoZero_NoOutliers_25Oct2023.csv", stringsAsFactors = T) %>% 
  dplyr::select(!(X))

###Zero read counts included
AllSpecies_AllMarkers_NoOutliers <- read.csv("AllSpecies_AllMarkers_NoOutliers_25Oct2023.csv", stringsAsFactors = T) %>% 
  dplyr::select(!(X))

###############Models with 0 read count removed###########
##############All Species################
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
hist(AllSpecies_AllMarkers_NoZero_NoOutliers$PropCorrectedReadsPerLitre) 
format(range(AllSpecies_AllMarkers_NoZero_NoOutliers$PropCorrectedReadsPerLitre),scientific=F) #[1] "0.000001240196" "0.999900000000"
var(AllSpecies_AllMarkers_NoZero_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.06758297
mean(AllSpecies_AllMarkers_NoZero_NoOutliers$PropCorrectedReadsPerLitre)# [1] 0.3001119


#GAM with cubic regression splines
#plot of data
Data_ScatterPlot <- ggplot(AllSpecies_AllMarkers_NoZero_NoOutliers,
                           aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre))+
  geom_point(aes(fill=Marker,shape=Marker),alpha=0.5)+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered")   + #label for the yaxis
  theme_bw()+
  theme(legend.position = "right",axis.text = element_text(size = 20), axis.title=element_text(size = 20),panel.border = element_rect(linewidth =0.5),plot.margin=unit(c(5,5,7,5), "mm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
Data_ScatterPlot

#get proper order of x values so the smooth function can be plotted properly
I1 <- order(AllSpecies_AllMarkers_NoZero_NoOutliers$QuantMeanPerLitre) 

#Make initial test GAM
Test1 <- gam(data = AllSpecies_AllMarkers_NoZero_NoOutliers,
            PropCorrectedReadsPerLitre ~ s(QuantMeanPerLitre, #Keep it simple with 1 variable for now
                                           fx=F,k=-1, # amount of smoothing not fixed to a preset value, cross-validation used to estimate optimal smoothing
                                           bs="cr")) #use cubic regression spline
summary(Test1)
anova(Test1)
vis.gam(Test1, theta=120, color="heat") #show the 3D version of the model
gam.check(Test1,pch=19,cex=.3)


#Plot Test1
plot(Test1, se=T) #basic plot of smoothing/GAM
#predict model values
Test1_pred <- predict(Test1, se=T, type="response")
#Plot for Test1 with data
Test1_ScatterPlot <- ggplot(AllSpecies_AllMarkers_NoZero_NoOutliers,
                            aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre))+
  geom_point(aes(fill=Marker,shape=Marker),alpha=0.5)+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  geom_line(aes(x=QuantMeanPerLitre[I1],  #main line
                y=Test1_pred$fit[I1]),
            colour="black")+
  geom_ribbon(aes(x = QuantMeanPerLitre[I1], 
                  ymin = Test1_pred$fit[I1]-2*Test1_pred$se[I1], #Lower CI
                  ymax = Test1_pred$fit[I1]+2*Test1_pred$se[I1]), #Upper CI
              fill = "grey", alpha = 0.5) +  
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="GAM: ReadsPerLitre ~ MeanCopiesPerLitre")   + #label for the yaxis
  theme_bw()+
  theme(legend.position = "right",axis.text = element_text(size = 20), axis.title=element_text(size = 20),panel.border = element_rect(linewidth =0.5),plot.margin=unit(c(5,5,7,5), "mm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
Test1_ScatterPlot

#make GAM with fill formula
Test2 <- gam(data = AllSpecies_AllMarkers_NoZero_NoOutliers,
             PropCorrectedReadsPerLitre ~ s(QuantMeanPerLitre,fx=F,k=-1,bs="cr") + Species)
               # Marker +
               # Type + 
               # Result + 
               # Species +
               # Run + 
               # Code+
summary(Test2)
anova(Test2)
vis.gam(Test2, theta=120, color="heat") #show the 3D version of the model
gam.check(Test2,pch=19,cex=.3)




#Plot Test2
plot(Test2, se=T) #basic plot of smoothing/GAM
#predict model values
Test2_pred <- predict(Test2, se=T, type="response")
#Plot for Test2 with data
Test2_ScatterPlot <- ggplot(AllSpecies_AllMarkers_NoZero_NoOutliers,
                            aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre))+
  geom_point(aes(fill=Marker,shape=Marker),alpha=0.5)+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
  scale_shape_manual(values=c(21,22,23))+
  geom_line(aes(x=QuantMeanPerLitre[I1],  #main line
                y=Test2_pred$fit[I1]),
            colour="black")+
  geom_ribbon(aes(x = QuantMeanPerLitre[I1], 
                  ymin = Test2_pred$fit[I1]-2*Test2_pred$se[I1], #Lower CI
                  ymax = Test2_pred$fit[I1]+2*Test2_pred$se[I1]), #Upper CI
              fill = "grey", alpha = 0.5) +  
  labs(x="Mean Copies Per Litre Filtered",
       y="Prop. Total Reads Per Litre Filtered",
       title="GAM: ReadsPerLitre ~ MeanCopiesPerLitre")   + #label for the yaxis
  theme_bw()+
  theme(legend.position = "right",axis.text = element_text(size = 20), axis.title=element_text(size = 20),panel.border = element_rect(linewidth =0.5),plot.margin=unit(c(5,5,7,5), "mm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
Test2_ScatterPlot







######save workspace######
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GAMs_workspace.RData")
