load("~/School/MASTER_Masters/R/Genetic_Environmental/LinearMixedModels_workspace.RData")

setwd("~/School/MASTER_Masters/Genetic_Oceano/LinearMixedModels_MultLinReg/")

library("nlme")
library("MuMIn")
library("lme4")
library("piecewiseSEM")
library("plyr")

#gen data --> use first PC of PCA on pop spec allele freqs?
#need a single response var

####Load data####
AllEnv_Outlier_FreqPCs <- read.csv("~/School/MASTER_Masters/Genetic_Data/June_2014_finals/PCA/PopFreqs/AllEnvNoMonth_Outlier_freqs_PCs.csv")
  rownames(AllEnv_Outlier_FreqPCs) <- AllEnv_Outlier_FreqPCs$Pop

CST_Outlier_FreqPCs <- read.csv("~/School/MASTER_Masters/Genetic_Data/June_2014_finals/PCA/PopFreqs/ChlASalTempNoMonth_Outlier_freqs_PCs.csv")
  rownames(CST_Outlier_FreqPCs) <- CST_Outlier_FreqPCs$Pop

QvalNorthNeutral <- read.csv("CorrelationNeutralQValue.csv")
  rownames(QvalNorthNeutral) <- QvalNorthNeutral$Pop
Dist_1D <- read.csv("~/School/MASTER_Masters/Sample_Sites/SampleSites_1DDistanceFromSUN.csv")
  rownames(Dist_1D) <- Dist_1D$Pop

AllEnv_EnvVars <- read.csv("~/School/MASTER_Masters/Oceanographic_Data/FINAL DATA/FINALOceano_Datasets/Oceano_data_stand_allvars_nomonths_trans_LMM_vars.csv")
  rownames(AllEnv_EnvVars) <- AllEnv_EnvVars$Pop

CST_EnvVars <- read.csv("~/School/MASTER_Masters/Oceanographic_Data/FINAL DATA/FINALOceano_Datasets/Oceano_data_stand_ChlATempSal_nomonths_trans_LMM_vars.csv")
  rownames(CST_EnvVars) <- CST_EnvVars$Pop


#all 12 temp vars
Temp_EnvVars <- read.csv("~/School/MASTER_Masters/Oceanographic_Data/FINAL DATA/FINALOceano_Datasets/Oceano_data_stand_allvars_nomonths_trans_AllTempLMM_vars.csv")
  rownames(Temp_EnvVars) <- Temp_EnvVars$Pop

#all 36 CST vars
CSTAll_EnvVars <- read.csv("~/School/MASTER_Masters/Oceanographic_Data/FINAL DATA/FINALOceano_Datasets/Oceano_data_stand_ChlATempSal_nomonths_trans.csv")
  rownames(CSTAll_EnvVars) <- CSTAll_EnvVars$Pop

#all 36 CST vars, not standardized
CST_noStand_EnvVars <- read.csv("~/School/MASTER_Masters/Oceanographic_Data/FINAL DATA/FINALOceano_Datasets/Oceano_data_ChlATempSal_nomonths_trans.csv")
  rownames(CST_noStand_EnvVars) <- CST_noStand_EnvVars$Pop

  

  
####Create DFs of gen*env for models####

#AllEnv, AllEnvOutPC1
AllEnv_LMM_data <- cbind(AllEnv_EnvVars, AllEnv_Outlier_FreqPCs$Axis1)
  colnames(AllEnv_LMM_data)[9] <- "AllEnvOutlierPCAxis1"
  AllEnv_LMM_data$Cluster <- c(rep("N", times=4), rep("S", times=8))
  AllEnv_LMM_data$Dist1D <- Dist_1D$X1DDist
  
#CST Env, CSTOutPC1
CST_LMM_data <- cbind(CST_EnvVars, CST_Outlier_FreqPCs$Axis1)
  colnames(CST_LMM_data)[14] <- "CSTOutlierPCAxis1"
  CST_LMM_data$Cluster <- c(rep("N", times=4), rep("S", times=8))
  CST_LMM_data$Dist1D <- Dist_1D$X1DDist
  
#Temp, AllEnvOut, CSTOut
AllEnvCST_Temp_LMM_data <- cbind(Temp_EnvVars, AllEnv_Outlier_FreqPCs$Axis1, CST_Outlier_FreqPCs$Axis1)
  colnames(AllEnvCST_Temp_LMM_data)[16:17] <- c("AllEnvOutlierPCAxis1","CSTOutlierPCAxis1")
  AllEnvCST_Temp_LMM_data$Cluster <- c(rep("N", times=4), rep("S", times=8))
  AllEnvCST_Temp_LMM_data$Dist1D <- Dist_1D$X1DDist

#AllEnvOut and CSTOut with all 36 CST vars
AllEnvCST_CSTAllEnvVars_LMs_data <- cbind(CSTAll_EnvVars, CST_Outlier_FreqPCs$Axis1, AllEnv_Outlier_FreqPCs$Axis1)
  colnames(AllEnvCST_CSTAllEnvVars_LMs_data)[40:41] <- c("CSTOutlierPCAxis1","AllEnvOutlierPCAxis1")
  AllEnvCST_CSTAllEnvVars_LMs_data$Cluster <- c(rep("N", times=4), rep("S", times=8))

#AllEnvOut and CSTOut with all 36 not standardized CST vars
AllEnvCST_CST_noStand_EnvVarss_LMs_data <- cbind(CST_noStand_EnvVars, CST_Outlier_FreqPCs$Axis1, AllEnv_Outlier_FreqPCs$Axis1)
  colnames(AllEnvCST_CST_noStand_EnvVarss_LMs_data)[40:41] <- c("CSTOutlierPCAxis1","AllEnvOutlierPCAxis1")
  AllEnvCST_CST_noStand_EnvVarss_LMs_data$Cluster <- c(rep("N", times=4), rep("S", times=8))

#####LMs on env vars vs Latitude, vs outlier axis####
SurfAvWinTemp <- lm(SurfAvWinTemp ~ Latitude, data = CST_all_EnvVars_LMs_data)
summary(lm(SurfAvWinTemp ~ Latitude, data = CST_all_EnvVars_LMs_data))$adj.r.squared

CST_allVars <- as.list(colnames(CST_all_EnvVars)[4:39])

Env_Lat_r2 <- list()
for (i in 4:39){
  EnvVar <- paste0(colnames(CST_all_EnvVars_LMs_data[i]))
  r2 <- summary(lm(CST_all_EnvVars_LMs_data[,i] ~ Latitude, data = CST_all_EnvVars_LMs_data))$adj.r.squared
  Env_Lat_r2[[EnvVar]] <- r2
  remove(EnvVar, r2)
}
Env_Lat_r2 <- as.data.frame(t(as.data.frame(Env_Lat_r2)))
Env_Lat_r2_Variable <- rownames(Env_Lat_r2)
Env_Lat_r2$Variable <- Env_Lat_r2_Variable
colnames(Env_Lat_r2)[1] <- "Latitude Adj R2"


Env_AllEnvOutPC1_r2 <- list()
for (i in 4:39){
  EnvVar <- paste0(colnames(CST_all_EnvVars_LMs_data[i]))
  r2 <- summary(lm(AllEnvOutlierPCAxis1 ~ CST_all_EnvVars_LMs_data[,i], data = CST_all_EnvVars_LMs_data))$adj.r.squared
  Env_AllEnvOutPC1_r2[[EnvVar]] <- r2
  remove(EnvVar, r2)
}
  Env_AllEnvOutPC1_r2 <- as.data.frame(t(as.data.frame(Env_AllEnvOutPC1_r2)))
  Env_AllEnvOutPC1_r2_Variable <- rownames(Env_AllEnvOutPC1_r2)
  Env_AllEnvOutPC1_r2$Variable <- Env_AllEnvOutPC1_r2_Variable
  colnames(Env_AllEnvOutPC1_r2)[1] <- "AllEnvOut Adj R2"


Env_CSTOutPC1_r2 <- list()
for (i in 4:39){
  EnvVar <- paste0(colnames(CST_all_EnvVars_LMs_data[i]))
  r2 <- summary(lm(CSTOutlierPCAxis1 ~ CST_all_EnvVars_LMs_data[,i], data = CST_all_EnvVars_LMs_data))$adj.r.squared
  Env_CSTOutPC1_r2[[EnvVar]] <- r2
  remove(EnvVar, r2)
}
  Env_CSTOutPC1_r2 <- as.data.frame(t(as.data.frame(Env_CSTOutPC1_r2)))
  Env_CSTOutPC1_r2_Variable <- rownames(Env_CSTOutPC1_r2)
  Env_CSTOutPC1_r2$Variable <- Env_CSTOutPC1_r2_Variable
  colnames(Env_CSTOutPC1_r2)[1] <- "CSTOut Adj R2"
  
  
Env_NoStand_AllEnvOutPC1_r2 <- list()
  for (i in 4:39){
    EnvVar <- paste0(colnames(CST_noStand_EnvVars_LMs_data[i]))
    r2 <- summary(lm(AllEnvOutlierPCAxis1 ~ CST_noStand_EnvVars_LMs_data[,i], data = CST_noStand_EnvVars_LMs_data))$adj.r.squared
    Env_NoStand_AllEnvOutPC1_r2[[EnvVar]] <- r2
    remove(EnvVar, r2)
  }
  Env_NoStand_AllEnvOutPC1_r2 <- as.data.frame(t(as.data.frame(Env_NoStand_AllEnvOutPC1_r2)))
  Env_NoStand_AllEnvOutPC1_r2_Variable <- rownames(Env_NoStand_AllEnvOutPC1_r2)
  Env_NoStand_AllEnvOutPC1_r2$Variable <- Env_NoStand_AllEnvOutPC1_r2_Variable
  colnames(Env_NoStand_AllEnvOutPC1_r2)[1] <- "NoStandAllEnvOut Adj R2"  
  
  
EnvNoStand_CSTOutPC1_r2 <- list()
  for (i in 4:39){
    EnvVar <- paste0(colnames(CST_noStand_EnvVars_LMs_data[i]))
    r2 <- summary(lm(CSTOutlierPCAxis1 ~ CST_noStand_EnvVars_LMs_data[,i], data = CST_noStand_EnvVars_LMs_data))$adj.r.squared
    EnvNoStand_CSTOutPC1_r2[[EnvVar]] <- r2
    remove(EnvVar, r2)
  }
EnvNoStand_CSTOutPC1_r2 <- as.data.frame(t(as.data.frame(EnvNoStand_CSTOutPC1_r2)))
  EnvNoStand_CSTOutPC1_r2_Variable <- rownames(EnvNoStand_CSTOutPC1_r2)
  EnvNoStand_CSTOutPC1_r2$Variable <- EnvNoStand_CSTOutPC1_r2_Variable
  colnames(EnvNoStand_CSTOutPC1_r2)[1] <- "NoStandCSTOut Adj R2"  
  
  

LM_AdjR2 <- join_all(list(Env_Lat_r2,Env_AllEnvOutPC1_r2,Env_CSTOutPC1_r2, 
                          EnvNoStand_CSTOutPC1_r2, Env_NoStand_AllEnvOutPC1_r2), 
                     by = 'Variable', type = 'full')
LM_AdjR2 <- LM_AdjR2[,c(2,1,3,4,5,6)]
LM_AdjR2 <- LM_AdjR2[1:36,]
View(LM_AdjR2)

write.csv(LM_AdjR2, "LM_AdjR2.csv")

#Pull out top 18 vars for each test (highest 18 R2) - adjusted in excel

LM_AdjR2 <- read.csv("LM_AdjR2.csv")

Latitude_list <- as.vector(LM_AdjR2[1:10,1])
AllEnvOutCST_list <- as.vector(LM_AdjR2[1:10,3])
CSTOutCST_list <- as.vector(LM_AdjR2[1:10,5])
AllEnvOutCSTNoStand_list <- as.vector(LM_AdjR2[1:10,7])
CSTOutCSTNoStand_list <- as.vector(LM_AdjR2[1:10,9])

Reduce(intersect, list(Latitude_list,AllEnvOutCST_list,CSTOutCST_list,AllEnvOutCSTNoStand_list,CSTOutCSTNoStand_list))  
#[1] "DepMinTemp"    "DepAvWinTemp"  "SurfMinTemp"   "SurfAvWinTemp" "DepAvSprTemp"  "DepAvAutSal"  

######Models#####
#use dredge for MuMIn?
#Response~Exp+Exp+Exp+Exp+(1|qvalue)
options(na.action = "na.fail")

AllEnv_GlobalModel <- lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp +
                            SurfAvAutSal + DepAvAutSal + DepMinSiO4)
  AllEnv_modelSel <- dredge(AllEnv_GlobalModel, evaluate=TRUE, trace=FALSE, beta="none")
  model.sel(AllEnv_modelSel) 

  #model averaging
  model.avg(AllEnv_modelSel, beta="none")

  

CST_GlobalModel <- glm(data =CST_LMM_data, CSTOutlierPCAxis1 ~ DepAvWinTemp + DepMinTemp + SurfAvWinTemp +
                         DepMaxSal + DepAvAutSal + SurfMaxChlA + SurfAvSprChlA + SurfAvSumChlA + DepMinChlA +
                         SurfMinChlA)
    CSTPC1
    DepAvWinTemp
    DepMinTemp
    SurfAvWinTemp
    DepMaxSal
    DepAvAutSal

    CSTPC4
    SurfMaxChlA
    SurfAvSprChlA
    SurfAvSumChlA
    DepMinChlA
    SurfMinChlA


  CST_modelSel <- dredge(CST_GlobalModel,evaluate = TRUE, trace=FALSE, beta="none")
  test <- model.sel(CST_modelSel)
  write.csv(as.data.frame(test),"CST_modelsel.csv") 
  
#   CST_model257 <- lme(data =CST_LMM_data, OutlierPCAxis1 ~ SurfAvWinTemp, random = ~1|QvalNorth)
#     AICc(CST_model257)
#     anova(CST_model257)
#     plot(CST_model257)
#     plot(CST_model257, OutlierPCAxis1 ~ SurfAvWinTemp)
#     summary(CST_model257)
#   
#   r.squaredGLMM(CST_model257)
  
    CST_modelNull <- glm(data =CST_LMM_data, CSTOutlierPCAxis1~1)
    AICc(CST_modelNull)

CST_HalfModel <- glm(data =CST_LMM_data, CSTOutlierPCAxis1 ~ DepAvWinTemp + DepMinTemp + SurfAvWinTemp +
                       DepMaxSal + DepAvAutSal)

  CSTHalf_modelSel <- dredge(CST_HalfModel,evaluate = TRUE, trace=FALSE, beta="none")
  model.sel(CSTHalf_modelSel)
  #model averaging
  model.avg(CSTHalf_modelSel, beta="none")
  

CST_BackHalfModel <- glm(data =CST_LMM_data, CSTOutlierPCAxis1 ~ SurfMaxChlA + SurfAvSprChlA + SurfAvSumChlA + DepMinChlA +
                           SurfMinChlA)

  CSTBackHalf_modelSel <- dredge(CST_BackHalfModel,evaluate = TRUE, trace=FALSE, beta="none")
  model.sel(CSTBackHalf_modelSel)
  model.avg(CSTBackHalf_modelSel, beta="none")
  
  
  
  
###issue - the global model has too many vars, model is being overdefined? or close to, can't
###pick out the actual important vars
###use just temp vars
AllEnv_Temp_GlobalModel <- glm(data =AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfAvSprTemp +SurfAvSumTemp +
                                 SurfAvAutTemp + DepAvWinTemp + DepAvSprTemp + DepAvSumTemp + DepAvAutTemp)

  AllEnv_Temp_modelSel <- dredge(AllEnv_Temp_GlobalModel, evaluate=TRUE, trace=FALSE, beta="none")
  model.sel(AllEnv_Temp_modelSel)
  #model averaging
  model.avg(AllEnv_Temp_modelSel, beta="none")
  

CST_Temp_GlobalModel <- glm(data =AllEnvCST_Temp_LMM_data, CSTOutlierPCAxis1 ~ SurfAvWinTemp + SurfAvSprTemp +SurfAvSumTemp +
                                   SurfAvAutTemp + DepAvWinTemp + DepAvSprTemp + DepAvSumTemp + DepAvAutTemp)

  CST_Temp_modelSel <- dredge(CST_Temp_GlobalModel, evaluate=TRUE, trace=FALSE, beta="none")
  model.sel(CST_Temp_modelSel)  
  #model averaging
  model.avg(CST_Temp_modelSel, beta="none")

  
  
    
  
  
  
  
  
  
  
  

  
  
  
  

####Models with most correlated vars####
#CST_all_EnvVars_LMs_data
# CSTOutlierPCAxis1 AllEnvOutlierPCAxis1
  
#use correlations from non-standardized data, but models with standardized?  
   
  
Reduce(intersect, list(Latitude_list,AllEnvOutCST_list,CSTOutCST_list,AllEnvOutCSTNoStand_list,CSTOutCSTNoStand_list))  
#[1] "DepMinTemp"    "DepAvWinTemp"  "SurfMinTemp"   "SurfAvWinTemp" "DepAvSprTemp"  "DepAvAutSal"  
Reduce(intersect, list(AllEnvOutCST_list,CSTOutCST_list,AllEnvOutCSTNoStand_list,CSTOutCSTNoStand_list))# - they all match!!!
# [1] "SurfAvWinTemp" "SurfMinTemp"   "SurfAvAutSal"  "DepMinChlA"    "SurfAvSumSal"  "DepMinTemp"    "DepAvWinTemp"  "SurfMinSal"    "DepAvAutSal"  
# [10] "DepAvSprTemp" 
Reduce(intersect, list(AllEnvOutCST_list[1:5],CSTOutCST_list[1:5],AllEnvOutCSTNoStand_list[1:5],CSTOutCSTNoStand_list[1:5]))
#[1] "SurfAvWinTemp" "SurfMinTemp"   "DepMinChlA"
Reduce(intersect, list(AllEnvOutCST_list[1:7],CSTOutCST_list[1:7],AllEnvOutCSTNoStand_list[1:7],CSTOutCSTNoStand_list[1:7])) 
#[1] "SurfAvWinTemp" "SurfMinTemp"   "SurfAvAutSal"  "DepMinChlA"    "DepMinTemp"    "DepAvWinTemp" 

#Response~Exp+Exp+Exp+Exp+(1|qvalue)
  
AllEnv_r2vars_GlobalModel <- lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp +
                                   SurfAvAutSal + DepMinChlA + DepMinTemp + DepAvWinTemp, 
                                 random = ~1|QvalNorth, method="ML")
  
  #since there are different fixed effects in each model, need REML to be FALSE (Zuur 2009)
  
  AllEnv_r2vars_modelSel <- dredge(AllEnv_r2vars_GlobalModel, evaluate=TRUE, trace=FALSE, beta="none")
  model.sel(AllEnv_r2vars_modelSel) 
  
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + SAS, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ SWT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ SMnT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + SWT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ SAS, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + SMnT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DWT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + SAS, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + DMnT + SAS, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DWT + DMnCA + SAS, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + SAS + SWT, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DWT + DMnCA, random = ~1|QvalNorth, method="ML"))
  r.squaredGLMM(lme(data =CST_all_EnvVars_LMs_data, AllEnvOutlierPCAxis1 ~ DMnCA + DMnT, random = ~1|QvalNorth, method="ML"))
  
  
CST_r2vars_GlobalModel <- lme(data =CST_all_EnvVars_LMs_data, CSTOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp +
                                SurfAvAutSal + DepMinChlA + DepMinTemp + DepAvWinTemp, 
                              random = ~1|QvalNorth, method="ML")

  CST_r2vars_modelSel <- dredge(CST_r2vars_GlobalModel, evaluate=TRUE, trace=FALSE, beta="none")
  model.sel(CST_r2vars_modelSel) 
  

#####Selected models and R2#####

  
  r.squared(glm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp +
        SurfAvAutSal + DepAvAutSal + DepMinSiO4))
  
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp + SurfAvAutSal + DepAvAutSal + DepMinSiO4))$adj.r.squared

#AllEnv PCA vars  
AllEnv_GlobalModel
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp + SurfAvAutSal + DepAvAutSal + DepMinSiO4))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutSal))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4 + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4 + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutSal + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4 + SurfAvAutSal))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutSal + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + SurfAvAutSal))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + DepMinSiO4))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4 + SurfAvWinTemp + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + SurfAvWinTemp + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4 + SurfAvAutSal + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + DepMinSiO4 + SurfAvAutSal))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + DepMinSiO4 + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepMinSiO4 + SurfAvAutSal + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutSal + SurfAvWinTemp + SurfMinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + SurfAvAutSal + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnv_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutSal + DepMinSiO4 + SurfMinTemp))$adj.r.squared


  
  
#AllEnv Out all temp
AllEnv_Temp_GlobalModel 
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~  SurfAvWinTemp + SurfAvSprTemp + SurfAvSumTemp + SurfAvAutTemp + DepAvWinTemp + DepAvSprTemp + DepAvSumTemp + DepAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvSprTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp + SurfAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + DepAvWinTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + DepAvWinTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + DepAvWinTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvSprTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + SurfAvAutTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + DepAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + DepAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + DepAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp + SurfAvSprTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvAutTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvAutTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + DepAvSumTemp))$adj.r.squared


  
#CSTOut all temp
CST_Temp_GlobalModel
summary(lm(data = AllEnvCST_Temp_LMM_data, CSTOutlierPCAxis1 ~ SurfAvWinTemp + SurfAvSprTemp + SurfAvSumTemp + SurfAvAutTemp + DepAvWinTemp + DepAvSprTemp + DepAvSumTemp + DepAvAutTem))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvSumTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvSprTemp + SurfAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + DepAvWinTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + DepAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvSumTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSprTemp + SurfAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + DepAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvSumTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvAutTemp + SurfAvAutTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ SurfAvSprTemp))$adj.r.squared
summary(lm(data = AllEnvCST_Temp_LMM_data, AllEnvOutlierPCAxis1 ~ DepAvWinTemp + SurfAvAutTemp + SurfAvWinTemp))$adj.r.squared

  
  

  
#####FINAL SAVE#####
save.image("~/School/MASTER_Masters/R/Genetic_Environmental/LinearMixedModels_workspace.RData")
