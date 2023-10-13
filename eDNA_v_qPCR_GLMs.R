#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")

#load libraries
library("tidyverse")
library("nlme")
library("MuMIn")
library("lme4")
library("stargazer")
library("ggeffects")


#set wd
setwd("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/")

#####load data (see eDNA_v_qPCR_DataProcessing.R)####
AtlSalmon_data <- read.csv("qPCR_eDNA_SalmoSalar_CleanedData_11Oct2023.csv") %>% 
  select(!(X))

PinkSalmon_data <-read.csv("qPCR_eDNA_PinkSalmon_CleanedData_11Oct2023.csv") %>% 
  select(!(X))

ArcticCharr_data <- read.csv("qPCR_eDNA_ArcticCharr_CleanedData_11Oct2023.csv") %>% 
  select(!(X)) 

####Data Prep####
#calculate QuantMean
AtlSalmon_data <- AtlSalmon_data %>% 
  mutate(QuantMean = rowMeans(select(.,Quant1,Quant2,Quant3), na.rm=TRUE)) %>% 
  relocate(QuantMean, .after=Quant3)

#Normalize CorrectedReads,CorrectedPercentReads,QuantMean by Vol filtered
AtlSalmon_data <- AtlSalmon_data %>% 
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

#####Data plots#####
# NormPercReads_12Steleo vs. NormQuantMean = P1
# NormCorrectedReads_12Steleo vs. NormQuantMean = P2
# log10(NormPercReads_12Steleo) vs. log10(NormQuantMean) = P3
# log10(NormCorrectedReads_12Steleo) vs. log10(NormQuantMean) = P4
# log10(NormPercReads_12Steleo) vs. log10(NormQuantMean) = P5
# log10(NormCorrectedReads_12Steleo) vs. log10(NormQuantMean) = P6


P1_base <- ggplot(AtlSalmon_data,aes(x=NormQuantMean,y=NormPercReads_12Steleo))
P1 <- P1_base +
  geom_point() +
  # scale_colour_manual(values = paletteer_d("ggthemes::Tableau_10"))+
  labs(x="Mean Copies\nPer Litre Filtered",
       y="% Total 12Steleo Reads\nPer Litre Filtered")   + #label for the yaxis
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P1 #changes the axes, etc

P2_base <- ggplot(AtlSalmon_data,aes(x=NormQuantMean,y=NormCorrectedReads_12Steleo))
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

P3_base <- ggplot(AtlSalmon_data,aes(x=log10(NormQuantMean),y=NormPercReads_12Steleo))
P3 <- P3_base +
  geom_point() +
  # scale_colour_manual(values = paletteer_d("ggthemes::Tableau_10"))+
    labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "% Total 12Steleo Reads\nPer Litre Filtered") +
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P3 #changes the axes, etc

P4_base <- ggplot(AtlSalmon_data,aes(x=log10(NormQuantMean),y=NormCorrectedReads_12Steleo))
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

P5_base <- ggplot(AtlSalmon_data,aes(x=log10(NormQuantMean),y=log10(NormPercReads_12Steleo)))
P5 <- P5_base +
  geom_point() +
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = expression(atop("log"["10"]~"(% Total 12Steleo Reads)", paste("Per Litre Filtered")))) +
  theme_bw() + #overall theme
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());P5 #changes the axes, etc

P6_base <- ggplot(AtlSalmon_data,aes(x=log10(NormQuantMean),y=log10(NormCorrectedReads_12Steleo)))
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
pdf("DataVisualization_Transformations_12Oct2023.pdf", height=20, width=15)
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
# First, separate out the variables I want for my models
AtlSalmon_model_data <- AtlSalmon_data %>% 
  select(SampleID,NormPercReads_12Steleo,NormQuantMean,log10NormQuantMean, RiverCode,Type,Year,Run,DNAConcScale,log10DNAConc)
#make dot charts of the continuous variables
op <- par(mfrow=c(3,2),mar=c(3,3,3,1))
dotchart(AtlSalmon_model_data$NormPercReads_12Steleo, main="NormPerReads",group=AtlSalmon_model_data$RiverCode)
plot(0,0,type="n", axes=F)
dotchart(AtlSalmon_model_data$NormQuantMean, main="MeanQuant",group=AtlSalmon_model_data$RiverCode)
dotchart(AtlSalmon_model_data$log10NormQuantMean, main="Log10MeanQuant",group=AtlSalmon_model_data$RiverCode)
dotchart(AtlSalmon_model_data$DNAConcScale, main="Scaled DNAConc",group=AtlSalmon_model_data$RiverCode)
dotchart(AtlSalmon_model_data$log10DNAConc, main="Log10DNAConc",group=AtlSalmon_model_data$RiverCode)
par(op)

#look at colinearity
z <- cbind(AtlSalmon_model_data$NormPercReads_12Steleo,
           AtlSalmon_model_data$NormQuantMean,
           AtlSalmon_model_data$log10NormQuantMean,
           AtlSalmon_model_data$DNAConcScale,
           AtlSalmon_model_data$log10DNAConc)
colnames(z) <- c("NormPerReads","NormQuantMean","log10Quant","DNAConcScale","Log10DNAConc")
pairs(z,
      lower.panel = panel.smooth,
      upper.panel=panel.cor,
      diag.panel=panel.hist)

corvif(z[,c(-1)])



#Vars
# River = Random
# Type WITHIN River = Random
# Year = Fixed # Random vars should have more than 5 levels
# Date same as year?
# Run = Random (control for effect of qPCR Run)
# DNA Conc = Random (Control for effect of DNA concentration) #Continuous var CANNOT be random, MUST be fixed

#full model
Full_LMM <- lmer(data = AtlSalmon_data, 
                 NormPercReads_12Steleo ~ log10NormQuantMean + DNAConc2 + Year + (1|Run) + (1|RiverCode/Type),
                 na.action = na.omit,
                 REML=T)
summary(Full_LMM)
AICc(Full_LMM)
stargazer(Full_LMM, 
          type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#sort out random using REML
Random1_LMM <- lmer(data = AtlSalmon_data, 
                    NormPercReads_12Steleo ~ log10NormQuantMean + DNAConc2 + Year + (1|RiverCode/Type),
                    na.action = na.omit,
                    REML=T)
Random2_LMM <- lmer(data = AtlSalmon_data, 
                    NormPercReads_12Steleo ~ log10NormQuantMean + DNAConc2 + Year + (1|Run),
                    na.action = na.omit,
                    REML=T)

#3 random models are: Full_LMM, Random1_LMM, Random2_LMM
summary(Full_LMM) #REML criterion at convergence: 2659.5
summary(Random1_LMM) #REML criterion at convergence: 2662.8
summary(Random2_LMM) #REML criterion at convergence: 2682.4

#Full model is best

#Select Fixed Effects, use ML and NOT REML for this step, as REML assumes fixed effects are correct
Full_LMM_ML <- lmer(data = AtlSalmon_data, 
                 NormPercReads_12Steleo ~ log10NormQuantMean + DNAConc2 + Year + (1|Run) + (1|RiverCode/Type),
                 na.action = na.omit,
                 REML=F)
Fixed1_LMM_ML <- lmer(data = AtlSalmon_data, 
                    NormPercReads_12Steleo ~ DNAConc2 + Year + (1|Run) + (1|RiverCode/Type),
                    na.action = na.omit,
                    REML=F)
Fixed2_LMM_ML <- lmer(data = AtlSalmon_data, 
                      NormPercReads_12Steleo ~ log10NormQuantMean + Year + (1|Run) + (1|RiverCode/Type),
                      na.action = na.omit,
                      REML=F)
Fixed3_LMM_ML <- lmer(data = AtlSalmon_data, 
                      NormPercReads_12Steleo ~ log10NormQuantMean + DNAConc2 + (1|Run) + (1|RiverCode/Type),
                      na.action = na.omit,
                      REML=F)
Fixed4_LMM_ML <- lmer(data = AtlSalmon_data, 
                    NormPercReads_12Steleo ~ log10NormQuantMean + (1|Run) + (1|RiverCode/Type),
                    na.action = na.omit,
                    REML=F)
Fixed5_LMM_ML <- lmer(data = AtlSalmon_data, 
                      NormPercReads_12Steleo ~ DNAConc2 + (1|Run) + (1|RiverCode/Type),
                      na.action = na.omit,
                      REML=F)
Fixed6_LMM_ML <- lmer(data = AtlSalmon_data, 
                      NormPercReads_12Steleo ~ Year + (1|Run) + (1|RiverCode/Type),
                      na.action = na.omit,
                      REML=F)

#7 fixed models are: Full_LMM_ML, Fixed1_LMM_ML, Fixed2_LMM_ML, Fixed3_LMM_ML, Fixed4_LMM_ML, Fixed5_LMM_ML, Fixed6_LMM_ML
AIC(Full_LMM_ML) #AIC 2698.728
AIC(Fixed1_LMM_ML) #AIC 2734.565
AIC(Fixed2_LMM_ML) #AIC 2697.943
AIC(Fixed3_LMM_ML) #AIC 2696.412
AIC(Fixed4_LMM_ML) #AIC 2695.762***
AIC(Fixed5_LMM_ML) #AIC 2733.943
AIC(Fixed6_LMM_ML) #AIC 2732.791
#Fixed4_LMM_ML is the best

#Final model
Final_LMM <- lmer(data = AtlSalmon_data, 
                      NormPercReads_12Steleo ~ log10NormQuantMean + (1|Run) + (1|RiverCode/Type),
                      na.action = na.omit,
                      REML=T)
summary(Final_LMM) #REML criterion at convergence: 2675.5
stargazer(Final_LMM, 
          type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


###plot the model
# Extract the prediction data frame
pred.mm <- ggpredict(Final_LMM, terms = c("log10NormQuantMean"))  # this gives overall predictions for the model

# Plot the predictions 
Final_LMM_plot <- ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data,                      # adding the raw data (scaled values)
             aes(x = log10(NormQuantMean), y = NormPercReads_12Steleo)) + 
  labs(x = expression(atop("log"["10"]~"(Mean Copies)", paste("Per Litre Filtered"))),
       y = "% Total 12Steleo Reads\nPer Litre Filtered") +
  theme_bw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank());Final_LMM_plot



 

#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
