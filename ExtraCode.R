#####################Extra##################
########## 12S Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

All_data_NoZero_Model_12S <- All_data_NoZero_Model %>% 
  filter(All_data_NoZero_Model$Marker == "12Steleo")

P1_base_12S_NoZero <- ggplot(All_data_NoZero_Model_12S,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_12S_NoZero <- P1_base_12S_NoZero +
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

P2_base_12S_NoZero <- ggplot(All_data_NoZero_Model_12S,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_12S_NoZero <- P2_base_12S_NoZero +
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

P3_base_12S_NoZero <- ggplot(All_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_12S_NoZero <- P3_base_12S_NoZero +
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

P4_base_12S_NoZero <- ggplot(All_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_12S_NoZero <- P4_base_12S_NoZero +
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

P5_base_12S_NoZero <- ggplot(All_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_12S_NoZero <- P5_base_12S_NoZero +
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

P6_base_12S_NoZero <- ggplot(All_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_12S_NoZero <- P6_base_12S_NoZero +
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
Multiplot_12S_NoZero <- grid.arrange(grobs=list(P1_12S_NoZero,P2_12S_NoZero,P3_12S_NoZero,P4_12S_NoZero,P5_12S_NoZero,P6_12S_NoZero),
                                     cols=2,
                                     top="Raw Data 12S_NoZero")
ggsave(Multiplot_12S_NoZero,
       file="DataVisualization_Transformations_12S_NoZero_20Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
All_data_NoZero_Prop_12S <- All_data_NoZero_Model_12S %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(All_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(All_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.0000009154" "1.0570512821"
boxplot(All_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre)
summary(All_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000009 0.0822997 0.2373737 0.2782800 0.4242000 1.0570513 
IQR <- IQR(All_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre)
quartiles <- quantile(All_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_Prop_12S_NoOutliers_Temp <- All_data_NoZero_Prop_12S %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_NoZero_Prop_12S_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(All_data_NoZero_Prop_12S_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#2 outliers removed

#explanatory variable - Mean Quant score
hist(All_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(All_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]   1.441469 248.628400
summary(All_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.441   7.564  12.400  17.991  20.885 248.628      28
boxplot(All_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(All_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(All_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

All_data_NoZero_Prop_12S_NoOutliers <- All_data_NoZero_Prop_12S_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(All_data_NoZero_Prop_12S_NoOutliers$QuantMeanPerLitre) 
boxplot(All_data_NoZero_Prop_12S_NoOutliers$QuantMeanPerLitre)
rm(All_data_NoZero_Prop_12S_NoOutliers_Temp)
#43 outliers removed


#redo plots with no outliers

P1_base_12S_NoZero_NoOut <- ggplot(All_data_NoZero_Prop_12S_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_12S_NoZero_NoOut <- P1_base_12S_NoZero_NoOut +
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

P3_base_12S_NoZero_NoOut <- ggplot(All_data_NoZero_Prop_12S_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_12S_NoZero_NoOut <- P3_base_12S_NoZero_NoOut +
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

P5_base_12S_NoZero_NoOut <- ggplot(All_data_NoZero_Prop_12S_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_12S_NoZero_NoOut <- P5_base_12S_NoZero_NoOut +
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
Multiplot_12S_NoZero_NoOut <- grid.arrange(grobs=list(P1_12S_NoZero_NoOut,P1_12S_NoZero,P2_12S_NoZero,
                                                      P3_12S_NoZero_NoOut,P3_12S_NoZero,P4_12S_NoZero,
                                                      P5_12S_NoZero_NoOut,P5_12S_NoZero,P6_12S_NoZero),
                                           cols=3,
                                           top="No Outliers 12S_NoZero")
ggsave(Multiplot_12S_NoZero_NoOut,
       file="DataVisualization_Transformations_12S_NoZero_NoOutliers_20Oct2023.pdf", 
       height=20, width=25,units = "in")





######### 12S_NoZero Final model ##########
FinalModel_12S_NoZero <- glmmTMB(data = charr_data_NoZero_Prop_12S_NoOutliers, 
                                 na.action = na.omit,family = beta_family(),REML=T,
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|Run) + (1|Code))
summary(FinalModel_12S_NoZero) 
AICc(FinalModel_12S_NoZero)# -314.5019
Anova(FinalModel_12S_NoZero)

###plot the model
# Extract the prediction data frame
pred.mm_12S_NoZero <- ggpredict(FinalModel_12S_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_12S_NoZero_plot <- ggplot(pred.mm_12S_NoZero) + 
  geom_line(data = pred.mm_12S_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_12S_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = charr_data_NoZero_Prop_12S_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM 12S_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_12S_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_12S_NoZero_GLMMBeta_20Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


########## FISHE Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

charr_data_NoZero_Model_FISHE <- charr_data_NoZero_Model %>% 
  filter(charr_data_NoZero_Model$Marker == "FISHE")

P1_base_FISHE_NoZero <- ggplot(charr_data_NoZero_Model_FISHE,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoZero <- P1_base_FISHE_NoZero +
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

P2_base_FISHE_NoZero <- ggplot(charr_data_NoZero_Model_FISHE,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_FISHE_NoZero <- P2_base_FISHE_NoZero +
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

P3_base_FISHE_NoZero <- ggplot(charr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoZero <- P3_base_FISHE_NoZero +
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

P4_base_FISHE_NoZero <- ggplot(charr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_FISHE_NoZero <- P4_base_FISHE_NoZero +
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

P5_base_FISHE_NoZero <- ggplot(charr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoZero <- P5_base_FISHE_NoZero +
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

P6_base_FISHE_NoZero <- ggplot(charr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_FISHE_NoZero <- P6_base_FISHE_NoZero +
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
Multiplot_FISHE_NoZero <- grid.arrange(grobs=list(P1_FISHE_NoZero,P2_FISHE_NoZero,P3_FISHE_NoZero,P4_FISHE_NoZero,P5_FISHE_NoZero,P6_FISHE_NoZero),
                                       cols=2,
                                       top="Raw Data FISHE_NoZero")
ggsave(Multiplot_FISHE_NoZero,
       file="DataVisualization_Transformations_FISHE_NoZero_20Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
charr_data_NoZero_Prop_FISHE <- charr_data_NoZero_Model_FISHE %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(charr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(charr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.005157143" "0.769381443"
boxplot(charr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
summary(charr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005157 0.096598 0.201000 0.260374 0.389711 0.769381 
IQR <- IQR(charr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
quartiles <- quantile(charr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

charr_data_NoZero_Prop_FISHE_NoOutliers_Temp <- charr_data_NoZero_Prop_FISHE %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#0 outliers removed

#explanatory variable - Mean Quant score
hist(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]     1.906399 83.533156
summary(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.906   8.117  13.422  16.838  21.222  83.533       3 
boxplot(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

charr_data_NoZero_Prop_FISHE_NoOutliers <- charr_data_NoZero_Prop_FISHE_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(charr_data_NoZero_Prop_FISHE_NoOutliers$QuantMeanPerLitre) 
boxplot(charr_data_NoZero_Prop_FISHE_NoOutliers$QuantMeanPerLitre)
rm(charr_data_NoZero_Prop_FISHE_NoOutliers_Temp)
#13 outliers removed


#redo plots with no outliers

P1_base_FISHE_NoZero_NoOut <- ggplot(charr_data_NoZero_Prop_FISHE_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoZero_NoOut <- P1_base_FISHE_NoZero_NoOut +
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

P3_base_FISHE_NoZero_NoOut <- ggplot(charr_data_NoZero_Prop_FISHE_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoZero_NoOut <- P3_base_FISHE_NoZero_NoOut +
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

P5_base_FISHE_NoZero_NoOut <- ggplot(charr_data_NoZero_Prop_FISHE_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoZero_NoOut <- P5_base_FISHE_NoZero_NoOut +
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
Multiplot_FISHE_NoZero_NoOut <- grid.arrange(grobs=list(P1_FISHE_NoZero_NoOut,P1_FISHE_NoZero,P2_FISHE_NoZero,
                                                        P3_FISHE_NoZero_NoOut,P3_FISHE_NoZero,P4_FISHE_NoZero,
                                                        P5_FISHE_NoZero_NoOut,P5_FISHE_NoZero,P6_FISHE_NoZero),
                                             cols=3,
                                             top="No Outliers FISHE_NoZero")
ggsave(Multiplot_FISHE_NoZero_NoOut,
       file="DataVisualization_Transformations_FISHE_NoZero_NoOutliers_20Oct2023.pdf", 
       height=20, width=25,units = "in")



######### FISHE_NoZero Final model ##########
FinalModel_FISHE_NoZero <- glmmTMB(data = charr_data_NoZero_Prop_FISHE_NoOutliers, 
                                   na.action = na.omit,family = beta_family(),REML=T,
                                   PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|Run) + (1|Code))
summary(FinalModel_FISHE_NoZero) 
AICc(FinalModel_FISHE_NoZero)# -185.9863
Anova(FinalModel_FISHE_NoZero)

###plot the model
# Extract the prediction data frame
pred.mm_FISHE_NoZero <- ggpredict(FinalModel_FISHE_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_FISHE_NoZero_plot <- ggplot(pred.mm_FISHE_NoZero) + 
  geom_line(data = pred.mm_FISHE_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_FISHE_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = charr_data_NoZero_Prop_FISHE_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM FISHE_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_FISHE_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_FISHE_NoZero_GLMMBeta_20Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width







########## MIFISHU Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

charr_data_NoZero_Model_MIFISHU <- charr_data_NoZero_Model %>% 
  filter(charr_data_NoZero_Model$Marker == "MIFISHU")

P1_base_MIFISHU_NoZero <- ggplot(charr_data_NoZero_Model_MIFISHU,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoZero <- P1_base_MIFISHU_NoZero +
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

P2_base_MIFISHU_NoZero <- ggplot(charr_data_NoZero_Model_MIFISHU,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_MIFISHU_NoZero <- P2_base_MIFISHU_NoZero +
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

P3_base_MIFISHU_NoZero <- ggplot(charr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoZero <- P3_base_MIFISHU_NoZero +
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

P4_base_MIFISHU_NoZero <- ggplot(charr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_MIFISHU_NoZero <- P4_base_MIFISHU_NoZero +
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

P5_base_MIFISHU_NoZero <- ggplot(charr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoZero <- P5_base_MIFISHU_NoZero +
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

P6_base_MIFISHU_NoZero <- ggplot(charr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_MIFISHU_NoZero <- P6_base_MIFISHU_NoZero +
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
Multiplot_MIFISHU_NoZero <- grid.arrange(grobs=list(P1_MIFISHU_NoZero,P2_MIFISHU_NoZero,P3_MIFISHU_NoZero,
                                                    P4_MIFISHU_NoZero,P5_MIFISHU_NoZero,P6_MIFISHU_NoZero),
                                         cols=2,
                                         top="Raw Data MIFISHU_NoZero")
ggsave(Multiplot_MIFISHU_NoZero,
       file="DataVisualization_Transformations_MIFISHU_NoZero_20Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
charr_data_NoZero_Prop_MIFISHU <- charr_data_NoZero_Model_MIFISHU %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(charr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(charr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.0000009191919" "0.9677906976744"
boxplot(charr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
summary(charr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000009 0.0621559 0.1945986 0.2559774 0.3913614 0.9677907 
IQR <- IQR(charr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
quartiles <- quantile(charr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp <- charr_data_NoZero_Prop_MIFISHU %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#7 outliers removed

#explanatory variable - Mean Quant score
hist(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]  1.138026 248.628400
summary(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.138   6.893  12.334  17.260  20.790 248.628      30 
boxplot(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

charr_data_NoZero_Prop_MIFISHU_NoOutliers <- charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(charr_data_NoZero_Prop_MIFISHU_NoOutliers$QuantMeanPerLitre) 
boxplot(charr_data_NoZero_Prop_MIFISHU_NoOutliers$QuantMeanPerLitre)
rm(charr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp)
#42 outliers removed


#redo plots with no outliers

P1_base_MIFISHU_NoZero_NoOut <- ggplot(charr_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoZero_NoOut <- P1_base_MIFISHU_NoZero_NoOut +
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

P3_base_MIFISHU_NoZero_NoOut <- ggplot(charr_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoZero_NoOut <- P3_base_MIFISHU_NoZero_NoOut +
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

P5_base_MIFISHU_NoZero_NoOut <- ggplot(charr_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoZero_NoOut <- P5_base_MIFISHU_NoZero_NoOut +
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
Multiplot_MIFISHU_NoZero_NoOut <- grid.arrange(grobs=list(P1_MIFISHU_NoZero_NoOut,P1_MIFISHU_NoZero,P2_MIFISHU_NoZero,
                                                          P3_MIFISHU_NoZero_NoOut,P3_MIFISHU_NoZero,P4_MIFISHU_NoZero,
                                                          P5_MIFISHU_NoZero_NoOut,P5_MIFISHU_NoZero,P6_MIFISHU_NoZero),
                                               cols=3,
                                               top="No Outliers MIFISHU_NoZero")
ggsave(Multiplot_MIFISHU_NoZero_NoOut,
       file="DataVisualization_Transformations_MIFISHU_NoZero_NoOutliers_20Oct2023.pdf", 
       height=20, width=25,units = "in")

######### MIFISHU_NoZero Final model ##########
FinalModel_MIFISHU_NoZero <- glmmTMB(data = charr_data_NoZero_Prop_MIFISHU_NoOutliers, 
                                     na.action = na.omit,family = beta_family(),REML=T,
                                     PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|Run) + (1|Code))

summary(FinalModel_MIFISHU_NoZero) 
AICc(FinalModel_MIFISHU_NoZero)# -340.6255
Anova(FinalModel_MIFISHU_NoZero)

###plot the model
# Extract the prediction data frame
pred.mm_MIFISHU_NoZero <- ggpredict(FinalModel_MIFISHU_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_MIFISHU_NoZero_plot <- ggplot(pred.mm_MIFISHU_NoZero) + 
  geom_line(data = pred.mm_MIFISHU_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_MIFISHU_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = charr_data_NoZero_Prop_MIFISHU_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM MIFISHU_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_MIFISHU_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_MIFISHU_NoZero_GLMMBeta_20Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width





###################Atlantic Salmon################
########## 12S Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

AtlSalmon_data_NoZero_Model_12S <- AtlSalmon_data_NoZero_Model %>% 
  filter(AtlSalmon_data_NoZero_Model$Marker == "12Steleo")

P1_base_12S_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_12S,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_12S_NoZero <- P1_base_12S_NoZero +
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

P2_base_12S_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_12S,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_12S_NoZero <- P2_base_12S_NoZero +
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

P3_base_12S_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_12S_NoZero <- P3_base_12S_NoZero +
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

P4_base_12S_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_12S_NoZero <- P4_base_12S_NoZero +
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

P5_base_12S_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_12S_NoZero <- P5_base_12S_NoZero +
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

P6_base_12S_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_12S,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_12S_NoZero <- P6_base_12S_NoZero +
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
Multiplot_12S_NoZero <- grid.arrange(grobs=list(P1_12S_NoZero,P2_12S_NoZero,P3_12S_NoZero,P4_12S_NoZero,P5_12S_NoZero,P6_12S_NoZero),
                                     cols=2,
                                     top="Raw Data 12S_NoZero")
ggsave(Multiplot_12S_NoZero,
       file="DataVisualization_Transformations_12S_NoZero_20Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_NoZero_Prop_12S <- AtlSalmon_data_NoZero_Model_12S %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.0000009154" "1.0570512821"
boxplot(AtlSalmon_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre)
summary(AtlSalmon_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000009 0.0822997 0.2373737 0.2782800 0.4242000 1.0570513 
IQR <- IQR(AtlSalmon_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_12S$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp <- AtlSalmon_data_NoZero_Prop_12S %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#2 outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]   1.441469 248.628400
summary(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.441   7.564  12.400  17.991  20.885 248.628      28
boxplot(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_12S_NoOutliers <- AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_12S_NoOutliers$QuantMeanPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_12S_NoOutliers$QuantMeanPerLitre)
rm(AtlSalmon_data_NoZero_Prop_12S_NoOutliers_Temp)
#43 outliers removed


#redo plots with no outliers

P1_base_12S_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_12S_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_12S_NoZero_NoOut <- P1_base_12S_NoZero_NoOut +
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

P3_base_12S_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_12S_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_12S_NoZero_NoOut <- P3_base_12S_NoZero_NoOut +
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

P5_base_12S_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_12S_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_12S_NoZero_NoOut <- P5_base_12S_NoZero_NoOut +
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
Multiplot_12S_NoZero_NoOut <- grid.arrange(grobs=list(P1_12S_NoZero_NoOut,P1_12S_NoZero,P2_12S_NoZero,
                                                      P3_12S_NoZero_NoOut,P3_12S_NoZero,P4_12S_NoZero,
                                                      P5_12S_NoZero_NoOut,P5_12S_NoZero,P6_12S_NoZero),
                                           cols=3,
                                           top="No Outliers 12S_NoZero")
ggsave(Multiplot_12S_NoZero_NoOut,
       file="DataVisualization_Transformations_12S_NoZero_NoOutliers_20Oct2023.pdf", 
       height=20, width=25,units = "in")





######### 12S_NoZero Final model ##########
FinalModel_12S_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_12S_NoOutliers, 
                                 na.action = na.omit,family = beta_family(),REML=T,
                                 PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|Run) + (1|Code))
summary(FinalModel_12S_NoZero) 
AICc(FinalModel_12S_NoZero)# -314.5019
Anova(FinalModel_12S_NoZero)

###plot the model
# Extract the prediction data frame
pred.mm_12S_NoZero <- ggpredict(FinalModel_12S_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_12S_NoZero_plot <- ggplot(pred.mm_12S_NoZero) + 
  geom_line(data = pred.mm_12S_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_12S_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_NoZero_Prop_12S_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#1B9E77"))+
  scale_shape_manual(values=c(21))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM 12S_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_12S_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_12S_NoZero_GLMMBeta_20Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


########## FISHE Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

AtlSalmon_data_NoZero_Model_FISHE <- AtlSalmon_data_NoZero_Model %>% 
  filter(AtlSalmon_data_NoZero_Model$Marker == "FISHE")

P1_base_FISHE_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_FISHE,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoZero <- P1_base_FISHE_NoZero +
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

P2_base_FISHE_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_FISHE,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_FISHE_NoZero <- P2_base_FISHE_NoZero +
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

P3_base_FISHE_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoZero <- P3_base_FISHE_NoZero +
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

P4_base_FISHE_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_FISHE_NoZero <- P4_base_FISHE_NoZero +
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

P5_base_FISHE_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoZero <- P5_base_FISHE_NoZero +
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

P6_base_FISHE_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_FISHE_NoZero <- P6_base_FISHE_NoZero +
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
Multiplot_FISHE_NoZero <- grid.arrange(grobs=list(P1_FISHE_NoZero,P2_FISHE_NoZero,P3_FISHE_NoZero,P4_FISHE_NoZero,P5_FISHE_NoZero,P6_FISHE_NoZero),
                                       cols=2,
                                       top="Raw Data FISHE_NoZero")
ggsave(Multiplot_FISHE_NoZero,
       file="DataVisualization_Transformations_FISHE_NoZero_20Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_NoZero_Prop_FISHE <- AtlSalmon_data_NoZero_Model_FISHE %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.005157143" "0.769381443"
boxplot(AtlSalmon_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
summary(AtlSalmon_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005157 0.096598 0.201000 0.260374 0.389711 0.769381 
IQR <- IQR(AtlSalmon_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp <- AtlSalmon_data_NoZero_Prop_FISHE %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#0 outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]     1.906399 83.533156
summary(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.906   8.117  13.422  16.838  21.222  83.533       3 
boxplot(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers <- AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers$QuantMeanPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers$QuantMeanPerLitre)
rm(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers_Temp)
#13 outliers removed


#redo plots with no outliers

P1_base_FISHE_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoZero_NoOut <- P1_base_FISHE_NoZero_NoOut +
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

P3_base_FISHE_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoZero_NoOut <- P3_base_FISHE_NoZero_NoOut +
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

P5_base_FISHE_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoZero_NoOut <- P5_base_FISHE_NoZero_NoOut +
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
Multiplot_FISHE_NoZero_NoOut <- grid.arrange(grobs=list(P1_FISHE_NoZero_NoOut,P1_FISHE_NoZero,P2_FISHE_NoZero,
                                                        P3_FISHE_NoZero_NoOut,P3_FISHE_NoZero,P4_FISHE_NoZero,
                                                        P5_FISHE_NoZero_NoOut,P5_FISHE_NoZero,P6_FISHE_NoZero),
                                             cols=3,
                                             top="No Outliers FISHE_NoZero")
ggsave(Multiplot_FISHE_NoZero_NoOut,
       file="DataVisualization_Transformations_FISHE_NoZero_NoOutliers_20Oct2023.pdf", 
       height=20, width=25,units = "in")



######### FISHE_NoZero Final model ##########
FinalModel_FISHE_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers, 
                                   na.action = na.omit,family = beta_family(),REML=T,
                                   PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|Run) + (1|Code))
summary(FinalModel_FISHE_NoZero) 
AICc(FinalModel_FISHE_NoZero)# -185.9863
Anova(FinalModel_FISHE_NoZero)

###plot the model
# Extract the prediction data frame
pred.mm_FISHE_NoZero <- ggpredict(FinalModel_FISHE_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_FISHE_NoZero_plot <- ggplot(pred.mm_FISHE_NoZero) + 
  geom_line(data = pred.mm_FISHE_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_FISHE_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_NoZero_Prop_FISHE_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM FISHE_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_FISHE_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_FISHE_NoZero_GLMMBeta_20Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width







########## MIFISHU Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

AtlSalmon_data_NoZero_Model_MIFISHU <- AtlSalmon_data_NoZero_Model %>% 
  filter(AtlSalmon_data_NoZero_Model$Marker == "MIFISHU")

P1_base_MIFISHU_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_MIFISHU,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoZero <- P1_base_MIFISHU_NoZero +
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

P2_base_MIFISHU_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_MIFISHU,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_MIFISHU_NoZero <- P2_base_MIFISHU_NoZero +
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

P3_base_MIFISHU_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoZero <- P3_base_MIFISHU_NoZero +
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

P4_base_MIFISHU_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_MIFISHU_NoZero <- P4_base_MIFISHU_NoZero +
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

P5_base_MIFISHU_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoZero <- P5_base_MIFISHU_NoZero +
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

P6_base_MIFISHU_NoZero <- ggplot(AtlSalmon_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_MIFISHU_NoZero <- P6_base_MIFISHU_NoZero +
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
Multiplot_MIFISHU_NoZero <- grid.arrange(grobs=list(P1_MIFISHU_NoZero,P2_MIFISHU_NoZero,P3_MIFISHU_NoZero,
                                                    P4_MIFISHU_NoZero,P5_MIFISHU_NoZero,P6_MIFISHU_NoZero),
                                         cols=2,
                                         top="Raw Data MIFISHU_NoZero")
ggsave(Multiplot_MIFISHU_NoZero,
       file="DataVisualization_Transformations_MIFISHU_NoZero_20Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
AtlSalmon_data_NoZero_Prop_MIFISHU <- AtlSalmon_data_NoZero_Model_MIFISHU %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(AtlSalmon_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(AtlSalmon_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.0000009191919" "0.9677906976744"
boxplot(AtlSalmon_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
summary(AtlSalmon_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000009 0.0621559 0.1945986 0.2559774 0.3913614 0.9677907 
IQR <- IQR(AtlSalmon_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp <- AtlSalmon_data_NoZero_Prop_MIFISHU %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#7 outliers removed

#explanatory variable - Mean Quant score
hist(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]  1.138026 248.628400
summary(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.138   6.893  12.334  17.260  20.790 248.628      30 
boxplot(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers <- AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers$QuantMeanPerLitre) 
boxplot(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers$QuantMeanPerLitre)
rm(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers_Temp)
#42 outliers removed


#redo plots with no outliers

P1_base_MIFISHU_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoZero_NoOut <- P1_base_MIFISHU_NoZero_NoOut +
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

P3_base_MIFISHU_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoZero_NoOut <- P3_base_MIFISHU_NoZero_NoOut +
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

P5_base_MIFISHU_NoZero_NoOut <- ggplot(AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoZero_NoOut <- P5_base_MIFISHU_NoZero_NoOut +
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
Multiplot_MIFISHU_NoZero_NoOut <- grid.arrange(grobs=list(P1_MIFISHU_NoZero_NoOut,P1_MIFISHU_NoZero,P2_MIFISHU_NoZero,
                                                          P3_MIFISHU_NoZero_NoOut,P3_MIFISHU_NoZero,P4_MIFISHU_NoZero,
                                                          P5_MIFISHU_NoZero_NoOut,P5_MIFISHU_NoZero,P6_MIFISHU_NoZero),
                                               cols=3,
                                               top="No Outliers MIFISHU_NoZero")
ggsave(Multiplot_MIFISHU_NoZero_NoOut,
       file="DataVisualization_Transformations_MIFISHU_NoZero_NoOutliers_20Oct2023.pdf", 
       height=20, width=25,units = "in")

######### MIFISHU_NoZero Final model ##########
FinalModel_MIFISHU_NoZero <- glmmTMB(data = AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers, 
                                     na.action = na.omit,family = beta_family(),REML=T,
                                     PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|Run) + (1|Code))

summary(FinalModel_MIFISHU_NoZero) 
AICc(FinalModel_MIFISHU_NoZero)# -340.6255
Anova(FinalModel_MIFISHU_NoZero)

###plot the model
# Extract the prediction data frame
pred.mm_MIFISHU_NoZero <- ggpredict(FinalModel_MIFISHU_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
FinalModel_MIFISHU_NoZero_plot <- ggplot(pred.mm_MIFISHU_NoZero) + 
  geom_line(data = pred.mm_MIFISHU_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_MIFISHU_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = AtlSalmon_data_NoZero_Prop_MIFISHU_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM MIFISHU_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(FinalModel_MIFISHU_NoZero_plot, #plot you want to save
       file = "AtlanticSalmon_MIFISHU_NoZero_GLMMBeta_20Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width





######Multiplot All FinalModel Plots#####
#Save all 6 plots in 1
Final_Multiplot_NoZero <- grid.arrange(grobs=list(FinalModel_NoZero_plot,
                                                  FinalModel_12S_NoZero_plot,
                                                  FinalModel_FISHE_NoZero_plot,
                                                  FinalModel_MIFISHU_NoZero_plot),
                                       cols=2,
                                       top="Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + (1|qPCRRun) + (1|River)")
ggsave(Final_Multiplot_NoZero,
       file="AtlanticSalmon_AllMarkers_NoZero_GLMMMultiplot_20Oct2023.pdf", 
       height=20, width=25,units = "in")




###################Arctic charr################
########## FISHE Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

ArcCharr_data_NoZero_Model_FISHE <- ArcCharr_data_NoZero_Model %>% 
  filter(ArcCharr_data_NoZero_Model$Marker == "FISHE")

P1_base_FISHE_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_FISHE,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoZero_charr <- P1_base_FISHE_NoZero_charr +
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

P2_base_FISHE_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_FISHE,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_FISHE_NoZero_charr <- P2_base_FISHE_NoZero_charr +
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

P3_base_FISHE_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoZero_charr <- P3_base_FISHE_NoZero_charr +
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

P4_base_FISHE_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_FISHE_NoZero_charr <- P4_base_FISHE_NoZero_charr +
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

P5_base_FISHE_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoZero_charr <- P5_base_FISHE_NoZero_charr +
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

P6_base_FISHE_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_FISHE,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_FISHE_NoZero_charr <- P6_base_FISHE_NoZero_charr +
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
Multiplot_FISHE_NoZero_charr <- grid.arrange(grobs=list(P1_FISHE_NoZero_charr,P2_FISHE_NoZero_charr,P3_FISHE_NoZero_charr,
                                                        P4_FISHE_NoZero_charr,P5_FISHE_NoZero_charr,P6_FISHE_NoZero_charr),
                                             cols=2,
                                             top="Raw Data FISHE_NoZero")
ggsave(Multiplot_FISHE_NoZero_charr,
       file="Charr_DataVisualization_Transformations_FISHE_NoZero_24Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
ArcCharr_data_NoZero_Prop_FISHE <- ArcCharr_data_NoZero_Model_FISHE %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(ArcCharr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(ArcCharr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.034050" "1.005051"
boxplot(ArcCharr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
summary(ArcCharr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03405 0.22458 0.65549 0.59030 0.94914 1.00505
IQR <- IQR(ArcCharr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre)
quartiles <- quantile(ArcCharr_data_NoZero_Prop_FISHE$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp <- ArcCharr_data_NoZero_Prop_FISHE %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#0 outliers removed

#explanatory variable - Mean Quant score
hist(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)# [1] 0.9816832 170.2106667
summary(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.9817   4.5617  10.8583  28.2332  31.7904 170.2107        2 
boxplot(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

ArcCharr_data_NoZero_Prop_FISHE_NoOutliers <- ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers$QuantMeanPerLitre) 
boxplot(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers$QuantMeanPerLitre)
rm(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers_Temp)
#11 outliers removed


#still have a proportion over 1, so remove that value
format(range(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers$PropCorrectedReadsPerLitre),scientific=F)
ArcCharr_data_NoZero_Prop_FISHE_NoOutliers <- ArcCharr_data_NoZero_Prop_FISHE_NoOutliers %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 1 row

#redo plots with no outliers

P1_base_FISHE_NoZero_NoOut_charr <- ggplot(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_FISHE_NoZero_NoOut_charr <- P1_base_FISHE_NoZero_NoOut_charr +
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

P3_base_FISHE_NoZero_NoOut_charr <- ggplot(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_FISHE_NoZero_NoOut_charr <- P3_base_FISHE_NoZero_NoOut_charr +
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

P5_base_FISHE_NoZero_NoOut_charr <- ggplot(ArcCharr_data_NoZero_Prop_FISHE_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_FISHE_NoZero_NoOut_charr <- P5_base_FISHE_NoZero_NoOut_charr +
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
Multiplot_FISHE_NoZero_NoOut_charr <- grid.arrange(grobs=list(P1_FISHE_NoZero_NoOut_charr,P1_FISHE_NoZero_charr,P2_FISHE_NoZero_charr,
                                                              P3_FISHE_NoZero_NoOut_charr,P3_FISHE_NoZero_charr,P4_FISHE_NoZero_charr,
                                                              P5_FISHE_NoZero_NoOut_charr,P5_FISHE_NoZero_charr,P6_FISHE_NoZero_charr),
                                                   cols=3,
                                                   top="No Outliers FISHE_NoZero")
ggsave(Multiplot_FISHE_NoZero_NoOut_charr,
       file="Charr_DataVisualization_Transformations_FISHE_NoZero_NoOutliers_24Oct2023.pdf", 
       height=20, width=25,units = "in")


#remove the plots once saved
rm(P1_charr,P1_base_charr,P1_base_NoOut_charr,P1_NoOut_charr,
   P2_charr,P2_base_charr,
   P3_charr,P3_base_charr,P3_base_NoOut_charr,P3_NoOut_charr,
   P4_charr,P4_base_charr,
   P5_charr,P5_base_charr,P5_base_NoOut_charr,P5_NoOut_charr,
   P6_charr,P6_base_charr)

######### FISHE_NoZero Final model ##########
Charr_FinalModel_FISHE_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_NoOutliers,
                                         PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|Code),
                                         na.action = na.omit,
                                         family = beta_family(),
                                         REML=T)
summary(Charr_FinalModel_FISHE_NoZero) 
AICc(Charr_FinalModel_FISHE_NoZero)# -151.7665
Anova(Charr_FinalModel_FISHE_NoZero, type="III")

###plot the model
# Extract the prediction data frame
pred.mm_Charr_FISHE_NoZero <- ggpredict(Charr_FinalModel_FISHE_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
Charr_FinalModel_FISHE_NoZero_plot <- ggplot(pred.mm_Charr_FISHE_NoZero) + 
  coord_cartesian(ylim=c(0, 1))+
  geom_line(data = pred.mm_Charr_FISHE_NoZero, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_Charr_FISHE_NoZero,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = ArcCharr_data_NoZero_Prop_FISHE_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#D95F02"))+
  scale_shape_manual(values=c(22))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM FISHE_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|River)") +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(Charr_FinalModel_FISHE_NoZero_plot, #plot you want to save
       file = "ArcticCharr_FISHE_NoZero_GLMMBeta_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width


########## MIFISHU Data plots and Outlier Removal##########
#picking data transformations
# PropCorrectedReadsPerLitre vs. QuantMeanPerLitre = P1
# CorrectedReadsPerLitre vs. QuantMeanPerLitre = P2
# PropCorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P3
# CorrectedReadsPerLitre vs. log10(QuantMeanPerLitre) = P4
# log10(PropCorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P5
# log10(CorrectedReadsPerLitre) vs. log10(QuantMeanPerLitre) = P6

ArcCharr_data_NoZero_Model_MIFISHU <- ArcCharr_data_NoZero_Model %>% 
  filter(ArcCharr_data_NoZero_Model$Marker == "MIFISHU")

P1_base_MIFISHU_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_MIFISHU,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoZero_charr <- P1_base_MIFISHU_NoZero_charr +
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

P2_base_MIFISHU_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_MIFISHU,aes(x=QuantMeanPerLitre,y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P2_MIFISHU_NoZero_charr <- P2_base_MIFISHU_NoZero_charr +
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

P3_base_MIFISHU_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoZero_charr <- P3_base_MIFISHU_NoZero_charr +
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

P4_base_MIFISHU_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=CorrectedReadsPerLitre,fill=Marker,shape=Marker))
P4_MIFISHU_NoZero_charr <- P4_base_MIFISHU_NoZero_charr +
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

P5_base_MIFISHU_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoZero_charr <- P5_base_MIFISHU_NoZero_charr +
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

P6_base_MIFISHU_NoZero_charr <- ggplot(ArcCharr_data_NoZero_Model_MIFISHU,aes(x=log10(QuantMeanPerLitre),y=log10(CorrectedReadsPerLitre),fill=Marker,shape=Marker))
P6_MIFISHU_NoZero_charr <- P6_base_MIFISHU_NoZero_charr +
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
Multiplot_MIFISHU_NoZero_charr <- grid.arrange(grobs=list(P1_MIFISHU_NoZero_charr,P2_MIFISHU_NoZero_charr,P3_MIFISHU_NoZero_charr,
                                                          P4_MIFISHU_NoZero_charr,P5_MIFISHU_NoZero_charr,P6_MIFISHU_NoZero_charr),
                                               cols=2,
                                               top="Raw Data MIFISHU_NoZero")
ggsave(Multiplot_MIFISHU_NoZero_charr,
       file="Charr_DataVisualization_Transformations_MIFISHU_NoZero_24Oct2023.pdf", 
       height=20, width=15,units = "in")



###remove outliers

#pull out just the Proportion data for modelling
ArcCharr_data_NoZero_Prop_MIFISHU <- ArcCharr_data_NoZero_Model_MIFISHU %>%  
  dplyr::select(!c(RawReads,RawReadsPerLitre,CorrectedReads,CorrectedReadsPerLitre))

#Response variable - proportion of total reads
hist(ArcCharr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre) #data is not normal, can't use z-score
format(range(ArcCharr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre),scientific=F)
# [1] "0.0000007895833" "1.0101010101010"
boxplot(ArcCharr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
summary(ArcCharr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.0000037 0.1412745 0.3462770 0.7124752 1.0101010 
IQR <- IQR(ArcCharr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre)
quartiles <- quantile(ArcCharr_data_NoZero_Prop_MIFISHU$PropCorrectedReadsPerLitre,probs=c(.25,.75),na.rm=F)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp <- ArcCharr_data_NoZero_Prop_MIFISHU %>% 
  filter(PropCorrectedReadsPerLitre > Lower & PropCorrectedReadsPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$PropCorrectedReadsPerLitre) 
boxplot(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$PropCorrectedReadsPerLitre)
#0 outliers removed

#explanatory variable - Mean Quant score
hist(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre) #data is not normal, can't use z-score
range(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
# [1]  0.6029412 170.2106667
summary(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.6029   3.1883   8.8512  26.3093  30.8081 170.2107       45
boxplot(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre)

IQR <- IQR(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,na.rm=T)
quartiles <- quantile(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp$QuantMeanPerLitre,probs=c(.25,.75),na.rm=T)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR

ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers <- ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp %>% 
  filter(QuantMeanPerLitre > Lower & QuantMeanPerLitre < Upper)
rm(IQR, quartiles,Lower,Upper)
hist(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers$QuantMeanPerLitre) 
boxplot(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers$QuantMeanPerLitre)
rm(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers_Temp)
#54 outliers removed

#still have a proportion over 1, so remove that value
format(range(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers$PropCorrectedReadsPerLitre),scientific=F)
ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers <- ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers %>% 
  filter(PropCorrectedReadsPerLitre < 1) #removed 2 rows

#redo plots with no outliers

P1_base_MIFISHU_NoZero_NoOut_charr <- ggplot(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=QuantMeanPerLitre,y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) #colour by marker type
P1_MIFISHU_NoZero_NoOut_charr <- P1_base_MIFISHU_NoZero_NoOut_charr +
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

P3_base_MIFISHU_NoZero_NoOut_charr <- ggplot(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=PropCorrectedReadsPerLitre,fill=Marker,shape=Marker))
P3_MIFISHU_NoZero_NoOut_charr <- P3_base_MIFISHU_NoZero_NoOut_charr +
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

P5_base_MIFISHU_NoZero_NoOut_charr <- ggplot(ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers,aes(x=log10(QuantMeanPerLitre),y=log10(PropCorrectedReadsPerLitre),fill=Marker,shape=Marker))
P5_MIFISHU_NoZero_NoOut_charr <- P5_base_MIFISHU_NoZero_NoOut_charr +
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
Multiplot_MIFISHU_NoZero_NoOut_charr <- grid.arrange(grobs=list(P1_MIFISHU_NoZero_NoOut_charr,P1_MIFISHU_NoZero_charr,P2_MIFISHU_NoZero_charr,
                                                                P3_MIFISHU_NoZero_NoOut_charr,P3_MIFISHU_NoZero_charr,P4_MIFISHU_NoZero_charr,
                                                                P5_MIFISHU_NoZero_NoOut_charr,P5_MIFISHU_NoZero_charr,P6_MIFISHU_NoZero_charr),
                                                     cols=3,
                                                     top="No Outliers MIFISHU_NoZero")
ggsave(Multiplot_MIFISHU_NoZero_NoOut_charr,
       file="Charr_DataVisualization_Transformations_MIFISHU_NoZero_NoOutliers_24Oct2023.pdf", 
       height=20, width=25,units = "in")

######### MIFISHU_NoZero Final model ##########
Charr_FinalModel_MIFISHU_NoZero <- glmmTMB(data = ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers, 
                                           PropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|Code),
                                           na.action = na.omit,
                                           family = beta_family(),
                                           REML=T)

summary(Charr_FinalModel_MIFISHU_NoZero) 
AICc(Charr_FinalModel_MIFISHU_NoZero)# -96.83043
Anova(Charr_FinalModel_MIFISHU_NoZero,type="III")

###plot the model
# Extract the prediction data frame
pred.mm_MIFISHU_NoZero_charr <- ggpredict(Charr_FinalModel_MIFISHU_NoZero, terms = c("QuantMeanPerLitre"))  # this gives overall predictions for the model

# Plot the predictions 
Charr_FinalModel_MIFISHU_NoZero_plot <- ggplot(pred.mm_MIFISHU_NoZero_charr) + 
  coord_cartesian(ylim=c(0, 1))+
  geom_line(data = pred.mm_MIFISHU_NoZero_charr, colour="black", 
            aes(x = x, y = predicted)) +          # slope
  geom_ribbon(data = pred.mm_MIFISHU_NoZero_charr,
              aes(x = x,
                  ymin = predicted - std.error,
                  ymax = predicted + std.error),
              fill = "grey", alpha = 0.5) +  # error band
  geom_point(data = ArcCharr_data_NoZero_Prop_MIFISHU_NoOutliers, alpha=0.7,         # adding the raw data (scaled values)
             aes(x = QuantMeanPerLitre, y = PropCorrectedReadsPerLitre,fill=Marker,shape=Marker)) +
  scale_fill_manual(values=c("#7570B3"))+
  scale_shape_manual(values=c(23))+
  labs(x = "Mean Copies Per Litre Filtered",
       y = "Prop. Total Reads Per Litre Filtered",
       title = "Beta GLMM MIFISHU_NoZero\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|River)") +
  # coord_cartesian(ylim=c(0, 1)) +
  theme_bw()+
  theme(legend.position = "right", 
        axis.text = element_text(size = 20), 
        axis.title=element_text(size = 20), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank())


#export plots
ggsave(Charr_FinalModel_MIFISHU_NoZero_plot, #plot you want to save
       file = "ArcticCharr_MIFISHU_NoZero_GLMMBeta_24Oct2023.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width





######Multiplot All FinalModel Plots#####
#Save all 6 plots in 1
Charr_Final_Multiplot_NoZero <- grid.arrange(grobs=list(Charr_FinalModel_NoZero_plot,
                                                        Charr_FinalModel_FISHE_NoZero_plot,
                                                        Charr_FinalModel_MIFISHU_NoZero_plot),
                                             cols=2,
                                             top="Beta GLMM\nPropCorrectedReadsPerLitre ~ QuantMeanPerLitre + DNAConcScale + (1|River)")
ggsave(Charr_Final_Multiplot_NoZero,
       file="ArcticCharr_AllMarkers_NoZero_GLMMMultiplot_24Oct2023.pdf", 
       height=20, width=10,units = "in")


###################Pink salmon################
