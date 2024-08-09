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

#t-test for two means of time points, 0 DHW and 15 DHW

t.test(TLC_CAM_log~Time_Point, data = rica)
t.test(TLC_CAM_log~ Time_Point, data = filter(rica, TLC_CAM_ratio < 2))

####################################################
## cor plot all data from AC DC 
####################################################
rica_acdc <- rica %>% select(Bleach_index, Bleach_index_rank, Heat_tolerance_CBASS,
                               symbiont_log, Growth_index, Growth_6mon, Source_Lat, 
                             TLC_CAM_log, delta_ratio, ddPCR_RICA_quiagen, Hcrit, Abundance) 

#plot the data frame
plot(rica_acdc)
cor(rica_acdc)

correlations <- cor(rica_ACDC, use="pairwise.complete.obs")

colnames(correlations) <- c( "Bleaching Index", 
                             "Bleaching Rank", "Heat Tolerance (CBASS)", "Symbiont Density",
                             "Specific Growth Rate", "6 Month Growth Rate", "Source Lat", "Aquarickettsia Adundance (ratio)",
                             "Delta Aquarickettsia (T1-T2)", "ddPCR_RICA_quiagen", "Hcrit")
rownames(correlations) <- c( "Bleaching Index", 
                             "Bleaching Rank", "Heat Tolerance (CBASS)", "Symbiont Density",
                             "Specific Growth Rate", "6 Month Growth Rate", "Source Lat", "Aquarickettsia Adundance (ratio)",
                             "Delta Aquarickettsia (T1-T2)", "ddPCR_RICA_quiagen", "Hcrit")
title1<- "AC DC Correlation Plot for T1 & T2"
cor_1<-corrplot(correlations,
         is.corr = F,
         method = "circle",
         col = COL1('YlOrRd', 200),
         title =title1,
         type = "upper",
         order = "original",
         tl.col = "black",
         tl.srt = 45,
         mar=c(0,0,1,0),
         addCoef.col = "black",
         tl.cex = 0.9)

####
#cor plot of t1
###
rica_ACDC2 <-  rica  %>% filter(Time_Point=="T1")%>% 
  select(Bleach_index, Bleach_index_rank, Heat_tolerance_CBASS,
         symbiont_log, Growth_index, Growth_6mon, Source_Lat, TLC_CAM_log, delta_ratio, ddPCR_RICA_quiagen, Hcrit) 

plot(rica_ACDC2)

correlations2 <- cor(rica_ACDC2, use="pairwise.complete.obs")

colnames(correlations2) <- c( "Bleaching Index", 
                             "Bleaching Rank", "Heat Tolerance (CBASS)", "Symbiont Density",
                             "Specific Growth Rate", "6 Month Growth Rate", "Source Lat", "Aquarickettsia Adundance (ratio)",
                             "Delta Aquarickettsia (T1-T2)", "ddPCR_RICA_quiagen", "Hcrit")
rownames(correlations2) <- c( "Bleaching Index", 
                             "Bleaching Rank", "Heat Tolerance (CBASS)", "Symbiont Density",
                             "Specific Growth Rate", "6 Month Growth Rate", "Source Lat", "Aquarickettsia Adundance (ratio)",
                             "Delta Aquarickettsia (T1-T2)", "ddPCR_RICA_quiagen", "Hcrit")

title2<- "AC DC Correlation Plot for T1"
cor_2<-corrplot(correlations2,
         is.corr = F,
         col = COL1('YlOrRd', 200),
         method = "circle",
         type = "upper",
         order = "original",
         tl.col = "black",
         tl.srt = 45,
         mar=c(0,0,1,0),
         addCoef.col = "black",
         tl.cex = 0.9,
         title =title2)


ggplot(rica)

####
#cor plot of t2
###
rica_ACDC3<-  rica  %>% filter(Time_Point=="T2")%>% 
  select(Bleach_index, Bleach_index_rank, Heat_tolerance_CBASS,
         symbiont_log, Growth_index, Growth_6mon, Source_Lat, TLC_CAM_log, delta_ratio, ddPCR_RICA_quiagen, Hcrit) 

correlations3 <- cor(rica_ACDC3, use="pairwise.complete.obs")

colnames(correlations3) <- c( "Bleaching Index", 
                              "Bleaching Rank", "Heat Tolerance (CBASS)", "Symbiont Density",
                              "Specific Growth Rate", "6 Month Growth Rate", "Source Lat", "Aquarickettsia Adundance (ratio)",
                              "Delta Aquarickettsia (T1-T2)", "ddPCR_RICA_quiagen", "Hcrit")
rownames(correlations3) <- c( "Bleaching Index", 
                              "Bleaching Rank", "Heat Tolerance (CBASS)", "Symbiont Density",
                              "Specific Growth Rate", "6 Month Growth Rate", "Source Lat", "Aquarickettsia Adundance (ratio)",
                              "Delta Aquarickettsia (T1-T2)", "ddPCR_RICA_quiagen", "Hcrit")

title3<- "AC DC Correlation Plot for T2"
cor_3<-corrplot(correlations3,
         is.corr = TRUE,
         col = COL1('YlOrRd', 200),
         method = "circle",
         type = "upper",
         order = "original",
         tl.col = "black",
         tl.srt = 45,
         mar=c(0,0,1,0),
         addCoef.col = "black",
         tl.cex = 0.9,
         title =title3)


#group all the cor plot together

ggarrange(cor_1, cor_3, cor_2, font.label = list(size = 10, color = "black"), labels=c("A", "B", "C"), common.legend = FALSE)

####################################
## Correlation tests for variables flagged in correlation plots
####################################

rica_t1 <- rica %>% filter(Time_Point=="T1")

#for rica levels and bleaching index
#0.2
cor.test(rica_t1$TLC_CAM_log, rica_t1$Bleach_index, 
        method ="pearson")

#for tlc abundance and growth index
cor.test(rica_t1$TLC_CAM_log, rica_t1$Growth_index, 
         method="pearson")

#for delta rica levels and Source_Lat
cor.test(rica_t1$TLC_CAM_log, rica_t1$Source_Lat, 
         method="pearson")

#for delta rica levels and Heat_tolerance
#0.3
cor.test(rica$delta_ratio, rica$Heat_tolerance_CBASS, 
         method="pearson")

#plot
ggplot(rica, aes(x=delta_ratio, y=Heat_tolerance_CBASS))+
  geom_point()


#for delta rica levels and symbiont_log
cor.test(rica_t1$delta_ratio, rica$symbiont_log, 
         method="pearson")


#for log ratio and dpcr data, just for fun and plot
cor.test(rica$TLC_CAM_log, rica$ddPCR_RICA_quiagen, 
         method="pearson")
#plot
ggplot(rica, aes(x=TLC_CAM_log, y=ddPCR_RICA_quiagen))+
  geom_point()

############################################
## Linear models with one predictor variable 
############################################

# growth model with RICA ratio being the response and growth rate being the predictor 
#Multiple R-squared:  0.03623
growthlm<- lm(TLC_CAM_log~Growth_index, data = rica)
summary(growthlm)

# bleaching index model with RICA ratio being the response and bleaching index being the predictor 
#Multiple R-squared:  0.002447
bleachlm<- lm(TLC_CAM_log~Bleach_index_rank, data = rica)
summary(bleachlm)

# symbiont density model with RICA ratio being the response and symbiont density  being the predictor 
#Multiple R-squared:  0.0001507
symbiontlm<- lm(TLC_CAM_log~symbiont_log, data = rica)
summary(symbiontlm)

# CBASS model with RICA ratio being the response and CBASS  being the predictor 
#Multiple R-squared:  0.002959
CBASSlm<-lm(TLC_CAM_log~Heat_tolerance_CBASS, data=rica)
summary(CBASSlm)

# Hcrit model with RICA ratio being the response and hcrit  being the predictor 
#Multiple R-squared:  0.01413
hcritlm<-lm(TLC_CAM_log~Hcrit, data=rica)
summary(hcritlm)

##################################################################
## Linear models with one predictor variable with the delta_ratio
##################################################################

rica_neg<- rica%>% filter(ratio_change=="neg")
rica_pos<- rica%>% filter(ratio_change=="pos")

# growth model with ratio being the response and growth rate being the predictor 
#Multiple R-squared:  0.2137
growthlm<- lm(TLC_CAM_log~Growth_index, data = rica_pos)
summary(growthlm)

#Multiple R-squared:  0.01341
growthlm2<- lm(TLC_CAM_log~Growth_index, data = rica_neg)
summary(growthlm2)

# bleaching index model with delta_ratio being the response and bleaching index being the predictor 
#Multiple R-squared:  0.1714,
bleachlm<- lm(TLC_CAM_log~Bleach_index_rank, data = rica_pos)
summary(bleachlm)

#Multiple R-squared:  0.01077
bleachlm2<- lm(TLC_CAM_log~Bleach_index_rank, data = rica_neg)
summary(bleachlm2)

# symbiont density model with delta_ratio being the response and symbiont density  being the predictor 
#Multiple R-squared:  0.1277
symbiontlm<- lm(TLC_CAM_log~symbiont_log, data = rica_pos)
summary(symbiontlm)

#Multiple R-squared:  0.002375
symbiontlm2<- lm(TLC_CAM_log~symbiont_log, data = rica_neg)
summary(symbiontlm2)

# CBASS model with delta_ratio being the response and CBASS  being the predictor 
#Multiple R-squared:  0.01411
CBASSlm<-lm(delta_ratio~Heat_tolerance_CBASS, data=rica)
summary(CBASSlm)

# Hcrit model with RICA ratio being the response and hcrit  being the predictor 
#Multiple R-squared:  0.01413
hcritlm<-lm(TLC_CAM_log~Hcrit, data=rica_pos)
summary(hcritlm)

#Multiple R-squared:  0.009735
hcritlm2<-lm(TLC_CAM_log~Hcrit, data=rica_neg)
summary(hcritlm2)

####
#using the rica abundance data as the predictor variables
####
##
##

rica_acdc <- rica %>%   select(Bleach_index, Bleach_index_rank, Heat_tolerance_CBASS,
         symbiont_log, Growth_index, Growth_6mon, Source_Lat, TLC_CAM_log, delta_ratio, ddPCR_RICA_quiagen, Hcrit, Abundance) 

library(MASS)
plot(rica_acdc)
#growth
full_G <- lm(Growth_index~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc)
stepG<-stepAIC(full_G, direction = "both", trace = FALSE)
summary(stepG)

glm_g <- glm(Growth_index~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc)
summary(glm_g)

anova(stepG)
glm_g$deviance


#Bleaching 
 
rica_acdc_b<- rica_acdc %>% filter(!is.na(Bleach_index)) %>%
  filter(!is.na(Abundance))%>%
  filter(!is.na(delta_ratio)) %>%
  filter(!is.na(ddPCR_RICA_quiagen))
  
full_B <- lm(Bleach_index~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance + delta_ratio * ddPCR_RICA_quiagen +
               ddPCR_RICA_quiagen * TLC_CAM_log + TLC_CAM_log * Abundance, data = rica_acdc_b)
stepB<-stepAIC(full_B, direction = "both", trace = FALSE)
summary(stepB)

glm_b <- glm(Bleach_index~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc)
summary(glm_b)

anova(stepB)
anova(glm_g)

#Heat_tolerance

rica_acdc_ht<- rica_acdc %>% filter(!is.na(Heat_tolerance_CBASS)) %>%
  filter(!is.na(Abundance))%>%
  filter(!is.na(delta_ratio)) %>%
  filter(!is.na(ddPCR_RICA_quiagen))

full_HT <- lm(Heat_tolerance_CBASS~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc_ht)
stepHT<-stepAIC(full_HT, direction = "both", trace = FALSE)
summary(stepHT) 

glm_ht <- glm(Heat_tolerance_CBASS~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc)
summary(glm_ht)

#symbiont density
rica_acdc_s<- rica_acdc %>% filter(!is.na(symbiont_log)) %>%
  filter(!is.na(Abundance))%>%
  filter(!is.na(delta_ratio)) %>%
  filter(!is.na(ddPCR_RICA_quiagen))

full_s <- lm(symbiont_log~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc_s)
steps<-stepAIC(full_s, direction = "both", trace = FALSE)
summary(steps)

glm_s <- glm(symbiont_log~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc)
summary(glm_s)

#hcrit 
rica_acdc_hc<- rica_acdc %>% filter(!is.na(Hcrit)) %>%
  filter(!is.na(Abundance))%>%
  filter(!is.na(delta_ratio)) %>%
  filter(!is.na(ddPCR_RICA_quiagen))

full_hc <- lm(Hcrit~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc_hc)
step_hc<-stepAIC(full_hc, direction = "both", trace = FALSE)
summary(step_hc)

glm_hc <- glm(Hcrit~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance, data = rica_acdc)
summary(glm_hc)

####################################################
## stepwise regression for all data from AC DC 
####################################################

full <- lm(TLC_CAM_log~., data = rica_ACDC)
null <- lm(TLC_CAM_log ~ 1, rica_ACDC)

#forward stepwise
forward<- step(null, direction = "forward", scope=formula(full), trace = 0)
forward$anova
summary(forward)


par(mfrow=c(2,2))
plot(forward)
dev.off()

#mass stepwise

library(MASS)

summary(stepAIC(null, direction = 'forward', scope = list(upper = full,
                                                          lower = null),trace = 0))


# different stepAIC
step<-stepAIC(full, direction = "both", trace = FALSE)

glm<-glm(formula = TLC_CAM_log ~ Bleach_index_rank + Growth_PropGrowth + 
           delta_ratio, data = rica_ACDC)
summary(step)
with(summary(step), 1 - deviance/null.deviance)

par(mfrow=c(2,2))
plot(step)
dev.off()


# no delta in the full model
rica_ACDC2<- rica_ACDC %>% select(TLC_CAM_log, Bleach_index_rank, Heat_tolerance_CBASS,
                                  symbiont_log, Growth_index) 

full <- lm(TLC_CAM_log~., data = rica_ACDC2)

step2<-stepAIC(full, direction = "both", trace = FALSE)

summary(step2)
with(summary(step2), 1 - deviance/null.deviance)

par(mfrow=c(2,2))
plot(step_hc)
dev.off()


#######
#testing the models
#######
set.seed(123)
train_samples<- rica_ACDC$TLC_CAM_log %>% createDataPartition(p=0.8, list = FALSE)
train<-rica_ACDC[train_samples,]
test<-rica_ACDC[-train_samples,]

pred<- step %>%predict(test)
data.frame(R2=R2(pred,test$TLC_CAM_log),
           RMSE=RMSE(pred,test$TLC_CAM_log),
           MAE=MAE(pred, test$TLC_CAM_log))

pred<- step2 %>%predict(test)
data.frame(R2=R2(pred,test$TLC_CAM_log),
           RMSE=RMSE(pred,test$TLC_CAM_log),
           MAE=MAE(pred, test$TLC_CAM_log))

pred<- growthlm %>%predict(test)
data.frame(R2=R2(pred,test$TLC_CAM_log),
           RMSE=RMSE(pred,test$TLC_CAM_log),
           MAE=MAE(pred, test$TLC_CAM_log))

#k fold

set.seed(123)
train.control<-trainControl(method = "repeatedcv", number = 10, repeats = 3)
stepK <- train(TLC_CAM_log~.,data=rica_ACDC, method="lm", trcontrol=train.control)
print(stepK)

train.control<-trainControl(method = "repeatedcv", number = 10, repeats = 3)
stepK2 <- train(TLC_CAM_log~.,data=rica_ACDC2, method="lm", trcontrol=train.control)
print(stepK2)

train.control<-trainControl(method = "repeatedcv", number = 10, repeats = 6)
stepK_grow<- train(TLC_CAM_log~Growth_PropGrowth,data=rica_ACDC2, method="lm", trcontrol=train.control)
print(stepK_grow)

###########################################
# model as delta as a predictor
###########################################

rica_ACDC_bi <-  rica  %>% select(ratio_change_bi, Bleach_index, Heat_tolerance_CBASS,
                                  symbiont_log, Growth_index) %>%
  filter(!is.na(Bleach_index)) %>%
  filter(!is.na(Heat_tolerance_CBASS)) %>%
  filter(!is.na(symbiont_log)) %>%
  filter(!is.na(Growth_index)) %>%
  filter(!is.na(ratio_change_bi))

bi_model<- glm(ratio_change_bi~
                 Bleach_index+
                 Heat_tolerance_CBASS+
                 Growth_index,
               family="binomial", data=rica_ACDC_bi)

summary(bi_model)

with(summary(bi_model), 1 - deviance/null.deviance)

par(mfrow=c(2,2))
plot(bi_model)

###
#ancova
###

delta_ <- lm(delta_ratio~., data = rica_ACDC)
anova(delta_)

#linear model to predict delta ratio

delta<- lm(delta_ratio~Heat_tolerance_CBASS+symbiont_log+Bleach_index, data = rica_ACDC)
summary(delta)

par(mfrow=c(2,2))
plot(delta)

#test the validity of the model
set.seed(123)
train.control<-trainControl(method = "repeatedcv", number = 10, repeats = 3)
deltak <- train(delta_ratio~Heat_tolerance_CBASS+symbiont_log+Bleach_index, data = rica_ACDC
                , method="lm", trcontrol=train.control)
print(deltak)


set.seed(123)
train_samples<- rica_ACDC$delta_ratio %>% createDataPartition(p=0.7, list = FALSE)
train<-rica_ACDC[train_samples,]
test<-rica_ACDC[-train_samples,]

pred<- delta %>%predict(test)
data.frame(R2=R2(pred,test$TLC_CAM_log),
           RMSE=RMSE(pred,test$TLC_CAM_log),
           MAE=MAE(pred, test$TLC_CAM_log))

####################################################
# Cluster Analysis
####################################################

##
#full cluster (all variables)
##
rica_cluster <-  rica  %>% select(Genotype, TLC_CAM_log, delta_ratio, Bleach_index_rank, 
                                  Heat_tolerance_CBASS, symbiont_log, Growth_index) %>%
  filter(!is.na(Bleach_index_rank)) %>%
  filter(!is.na(Heat_tolerance_CBASS)) %>%
  filter(!is.na(symbiont_log)) %>%
  filter(!is.na(Growth_PropGrowth)) %>%
  filter(!is.na(delta_ratio))

#label the rownames with genotype

rica_cluster1<- rica_cluster[,-1]
rica_matrix <- as.matrix(rica_cluster1)
rownames(rica_matrix) <- rica_cluster[,1]
rica_cluster<- as.data.frame(rica_matrix)
head(rica_cluster)

#average cluster
rica_distance <- dist(rica_cluster, method = "euclidean")
rica_hc_single <- hclust(rica_distance, method= "average")
fviz_dend(rica_hc_single, cex = .5, k=7, palette = NULL)

##
#just growth, ratio and change
##

rica_cluster2<- rica %>% select(Genotype, Growth_PropGrowth, delta_ratio, TLC_CAM_log)%>% 
  filter(!is.na(Growth_PropGrowth))%>%
  filter(!is.na(delta_ratio))

#label the rownames with genotype
rica_cluster3<- rica_cluster2[,-1]
rica_matrix2 <- as.matrix(rica_cluster3)
rownames(rica_matrix2) <- rica_cluster2[,1]
rica_cluster3<- as.data.frame(rica_matrix2)
head(rica_cluster3)

#average cluster 
rica_distance2 <- dist(rica_cluster3, method = "euclidean")
rica_hc_single2 <- hclust(rica_distance2, method= "average")
fviz_dend(rica_hc_single2, cex = .5, k=5, palette = NULL)

#####
## just ratio and delta
####

rica_cluster<- rica %>% filter(DHW==18) %>% select(Genotype, TLC_CAM_log)

#label the rownames with genotype
rica_cluster1<- rica_cluster[,-1]
rica_matrix <- as.matrix(rica_cluster1)
rownames(rica_matrix) <- rica_cluster[,1]
rica_cluster<- as.data.frame(rica_matrix)
head(rica_cluster)

#average cluster 
rica_distance <- dist(rica_cluster, method = "euclidean")
rica_hc_single <- hclust(rica_distance, method= "average")
fviz_dend(rica_hc_single, cex = .5, k=5, palette = NULL)


#######################################################
###cluster heat maps
#######################################################

rica_ACDC_scaled <- scale(rica_cluster)
pheatmap(t(rica_ACDC_scaled), cutree_cols=4)

rica_ACDC_scaled2 <- scale(rica_cluster3)
pheatmap(t(rica_ACDC_scaled2), cutree_cols=10)

#######################################################
## PCA
########################################################

rica_cov <-cov(rica_ACDC)

pca <- prcomp(rica_ACDC, scale. = FALSE)
summary(pca)
head(rica_ACDC)

fviz_pca_biplot(pca, repel=TRUE)

#just growth, delta, and ratio
rica3<- rica %>% filter(!is.na(Growth_index))%>%
  filter(TLC_CAM_ratio < 2) %>%
  filter(ratio_change !="NA")

rica_pca<- rica3 %>% select(Growth_index, TLC_CAM_log, delta_ratio)%>%
  filter(!is.na(delta_ratio))

rica_label1<- rica3 %>% select(Growth_index, TLC_CAM_log, Genotype, Time_Point, delta_ratio)%>%
  filter(!is.na(delta_ratio))

pca2 <- prcomp(rica_pca, scale. = FALSE)
summary(pca2)
head(rica_ACDC)

fviz_pca(pca2, axes = c(1, 2), habillage=rica_label1$Time_Point,   repel=TRUE)


###
#SEM models 
###


fit1 <- brm(formula = symbiont_log~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance,
            data = rica)

summary(fit1)

#plot the model
pp_check(fit1)
#plot the model
bayesplot::ppc_scatter_avg(y = rica_acdc_s$symbiont_log, yrep = posterior_predict(fit1))


#make a plot with the model
post_sampels <- as.data.frame(fit1)

ggplot(post_sampels, aes(x=b_delta_ratio))+ 
  geom_histogram(bins = 30, fill="blue")


##
fit2 <- brm(formula = Bleach_index~ delta_ratio + ddPCR_RICA_quiagen + TLC_CAM_log + Abundance,
            data = rica)

summary(fit2)

#plot the model
pp_check(fit2)


#plot the model
bayesplot::ppc_scatter_avg(y = rica_acdc_b$Bleach_index, yrep = posterior_predict(fit2))

#make a plot with the model
post_sampels <- as.data.frame(fit2)

ggplot(post_sampels, aes(x=b_Abundance))+ 
  geom_histogram(bins = 30, fill="blue")



