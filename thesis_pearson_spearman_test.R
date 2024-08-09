library(tidyverse) 
library(devtools)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(forcats)
library(pheatmap)
library(car)
library(caret)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(brms)
library(bayesplot)


#######
## Sterling Butler 
## NOAA, RSMAS University of Miami 
## 12/14/2023
## Aquarickettsia x Acropora cervicornis  
## UM nursery models and statistical anaylsis
#######

#set working directory
setwd("/Volumes/RSMASFILES2/RICA_DATA_NEW/")

#get metadata
rica <- read.csv("06262024_RICA_Nursery_qPCR_META.csv")
head(rica)

### filter out outliers
rica<- rica %>% filter(TLC_CAM_ratio < 2)

#add the groups
rica<- rica %>%
  mutate(RICA_Abundance = case_when(TLC_CAM_ratio < 0.7652562 ~ 'LOW',
                                    TLC_CAM_ratio < 1.3418562 ~ 'MED',
                                    TLC_CAM_ratio < 2 ~ 'HIGH'))

#test for normalility

shapiro.test(rica$TLC_CAM_log)
shapiro.test(rica$delta_ratio)
shapiro.test(rica$symbiont_log)
shapiro.test(rica$Bleach_index)
shapiro.test(rica$Growth_index)
shapiro.test(rica$Growth_Specific)
shapiro.test(rica$Heat_tolerance_CBASS)

###Look into transformation for this one
shapiro.test(rica$Growth_6mon)
hist(rica$Growth_6mon)

####################################
## Correlation tests for t1 table 
####################################

rica_t1 <- rica %>% filter(Time_Point=="T1")

#for rica levels and bleaching index
cor.test(rica_t1$TLC_CAM_log, rica_t1$Bleach_index, 
         method ="pearson")
cor.test(rica_t1$TLC_CAM_log, rica_t1$Bleach_index, 
          method ="spearman")

#for rica levels and heat tolerance index
cor.test(rica_t1$TLC_CAM_log, rica_t1$Heat_tolerance_CBASS, 
         method ="pearson")
cor.test(rica_t1$TLC_CAM_log, rica_t1$Heat_tolerance_CBASS, 
         method ="spearman")

#for rica levels and symbiont_log
cor.test(rica_t1$TLC_CAM_log, rica_t1$symbiont_log, 
         method ="pearson")
cor.test(rica_t1$TLC_CAM_log, rica_t1$symbiont_log, 
         method ="spearman")

#for rica levels and Growth_6mon
cor.test(rica_t1$TLC_CAM_log, rica_t1$Growth_6mon, 
         method ="pearson")
cor.test(rica_t1$TLC_CAM_log, rica_t1$Growth_6mon, 
         method ="spearman")

#for rica levels and Growth_index
cor.test(rica_t1$TLC_CAM_log, rica_t1$Growth_index, 
         method ="pearson")
cor.test(rica_t1$TLC_CAM_log, rica_t1$Growth_index, 
         method ="spearman")

#for rica levels and Growth_index
cor.test(rica_t1$TLC_CAM_log, rica_t1$Source_Lat, 
         method ="pearson")
cor.test(rica_t1$TLC_CAM_log, rica_t1$Source_Lat, 
         method ="spearman")

####################################
## Correlation tests for t1 + T2 table 
####################################



#for rica levels and bleaching index
cor.test(rica$TLC_CAM_log, rica$Bleach_index, 
         method ="pearson")
cor.test(rica$TLC_CAM_log, rica$Bleach_index, 
         method ="spearman")

#for rica levels and heat tolerance index
cor.test(rica$TLC_CAM_log, rica$Heat_tolerance_CBASS, 
         method ="pearson")
cor.test(rica$TLC_CAM_log, rica$Heat_tolerance_CBASS, 
         method ="spearman")

#for rica levels and symbiont_log
cor.test(rica$TLC_CAM_log, rica$symbiont_log, 
         method ="pearson")
cor.test(rica$TLC_CAM_log, rica$symbiont_log, 
         method ="spearman")

#for rica levels and Growth_6mon
cor.test(rica$TLC_CAM_log, rica$Growth_6mon, 
         method ="pearson")
cor.test(rica$TLC_CAM_log, rica$Growth_6mon, 
         method ="spearman")

#for rica levels and Growth_index
cor.test(rica$TLC_CAM_log, rica$Growth_index, 
         method ="pearson")
cor.test(rica$TLC_CAM_log, rica$Growth_index, 
         method ="spearman")

#for rica levels and Growth_index
cor.test(rica$TLC_CAM_log, rica$Source_Lat, 
         method ="pearson")
cor.test(rica$TLC_CAM_log, rica$Source_Lat, 
         method ="spearman")

####################################
## Correlation tests for delta ratio
####################################

#for rica levels and bleaching index
cor.test(rica$delta_ratio, rica$Bleach_index, 
         method ="pearson")
cor.test(rica$delta_ratio, rica$Bleach_index, 
         method ="spearman")

#for rica levels and heat tolerance index
cor.test(rica$delta_ratio, rica$Heat_tolerance_CBASS, 
         method ="pearson")
cor.test(rica$delta_ratio, rica$Heat_tolerance_CBASS, 
         method ="spearman")

#for rica levels and symbiont_log
cor.test(rica$delta_ratio, rica$symbiont_log, 
         method ="pearson")
cor.test(rica$delta_ratio, rica$symbiont_log, 
         method ="spearman")

#for rica levels and Growth_6mon
cor.test(rica$delta_ratio, rica$Growth_6mon, 
         method ="pearson")
cor.test(rica$delta_ratio, rica$Growth_6mon, 
         method ="spearman")

#for rica levels and Growth_index
cor.test(rica$delta_ratio, rica$Growth_index, 
         method ="pearson")
cor.test(rica$delta_ratio, rica$Growth_index, 
         method ="spearman")

#for rica levels and Growth_index
cor.test(rica$delta_ratio, rica$Source_Lat, 
         method ="pearson")
cor.test(rica$delta_ratio, rica$Source_Lat, 
         method ="spearman")





