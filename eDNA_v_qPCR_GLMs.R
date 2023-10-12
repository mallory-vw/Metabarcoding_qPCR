#load workspace
load("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")

#load libraries
library("tidyverse")
library("nlme")
library("MuMIn")
library("lme4")
library("merTools")
library("colorBlindness")


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
  relocate(NormQuantMean, .after=QuantMean)

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
testlm <- lm(NormPercReads_12Steleo ~ log10(NormQuantMean),
              data = AtlSalmon_data)

anova(testlm)
summary(testlm)








options(na.action = "na.fail")

AllEnv_GlobalModel <- lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp +
                           SurfAvAutSal + DepAvAutSal + DepMinSiO4)
AllEnv_modelSel <- dredge(AllEnv_GlobalModel, evaluate=TRUE, trace=FALSE, beta="none")
model.sel(AllEnv_modelSel) 

#model averaging
model.avg(AllEnv_modelSel, beta="none")






#model variable selections

 




#save workspace
save.image("C:/Users/vanwyngaardenma/Documents/Bradbury/Metabarcoding/Metabardoding_qPCR/eDNA_v_qPCR_GLMs_workspace.RData")
