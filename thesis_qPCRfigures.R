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
library(ggpubr)

#######
## Sterling Butler 
## NOAA, RSMAS University of Miami 
## 12/04/2023
## Aquarickettsia x Acropora cervicornis  
## UM nursery data manipulation and exploratory graphs
#######

#set working directory
setwd("/Volumes/RSMASFILES2/RICA_DATA_NEW/")

#get metadata
rica <- read.csv("08022024_RICA_Nursery_qPCR_META.csv")

head(rica)

### adding a column for positive or negative delta ratio changes
rica <-rica %>%
  mutate(ratio_change = if_else(delta_ratio > 0, "+", "-"))

#add the groups
rica<- rica %>%
  mutate(RICA_Abundance = case_when(TLC_CAM_ratio < 0.7652562 ~ 'LOW',
                                    TLC_CAM_ratio < 1.3418562 ~ 'MED',
                                    TLC_CAM_ratio < 2 ~ 'HIGH'))

head(rica)

### filter out outlier
rica<- rica %>% filter(ratio_mean != "NA")

#if you want to show what is significant and then what you want to compare
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
comparisons_material = list( c("T1", "T2"))

qpcr_drop<- ggplot(rica, aes(x=Time_Point, y=TLC_CAM_log))+
  geom_boxplot(alpha = 2, aes(colour=Time_Point))+
  theme_classic2()+
  scale_color_manual(values = c("#50a0aa", "#aa5a50"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  labs(title = "qPCR", 
       x = "Time Point", y = "log Tlc1 to CaM Ratio ",
       color= "Time Point")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)


ggsave("/Volumes/RSMASFILES2/ThesisPlots/fig_1.png", width = 5.5, height = 4.5, units = "in",
       dpi=300)

summary(aov(TLC_CAM_log~Time_Point, rica))


drop_16s <-ggplot(rica, aes(x=Time_Point, y=rica_16S_abund))+
  geom_boxplot(alpha = 2, aes(colour=Time_Point))+
  theme_classic2()+
  scale_color_manual(values = c("#50a0aa", "#aa5a50"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  labs(title = "16S rRNA", 
       x = "Time Point", y = "MD3-55 Relative Abundance",
       color= "Time Point")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)


drop_plot <- ggarrange(qpcr_drop, drop_16s, font.label = list(size = 10, color = "black"), labels=c("A", "B"), common.legend = TRUE)

ggsave("/Volumes/RSMASFILES2/ThesisPlots/abundanceBox.png", width = 8, height = 6, units = "in",
       dpi=300)

wil_16 <-wilcox.test(rica_16S_abund ~ Time_Point, data = rica,  exact = FALSE)

wil_qpcr <-wilcox.test(TLC_CAM_log ~ Time_Point, data = rica,  exact = FALSE)

#mean and sd for timepoints
#t1
rica_t1<- rica %>% filter(Time_Point == "T1")
t1_mean<- round(mean(rica_t1$TLC_CAM_ratio), 3)
t1_sd<- round(sd(rica_t1$TLC_CAM_ratio), 3)
#t2
rica_t2<- rica %>% filter(Time_Point == "T2")
t2_mean<- round(mean(rica_t2$TLC_CAM_ratio), 3)
t2_sd<- round(sd(rica_t2$TLC_CAM_ratio), 3)

# Boxplot of RICA levels for time_points, separated on the x axis by
# Genotypex=fct_infreq(Time_Point)

facet_g <- ggplot(rica, aes(x=Time_Point, y=TLC_CAM_log))+
  geom_boxplot(aes(group=Genotype_Timepoint, colour=Time_Point))+
  facet_wrap(~Genotype_a)+
  theme_classic2()+
  scale_color_manual(values = c("#50a0aa", "#aa5a50"))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 6),
        title = element_text(size = 12))+
  labs( x = "Time Point", y = "log tlc1 to CaM ratio",  color="Time Point") 


#####
## full boxplot for all genotypes, with color for ratio change
#####

full_g <- ggplot(rica, aes(x=fct_reorder(Genotype_a, TLC_CAM_log), y=TLC_CAM_log, color= ratio_change))+
  geom_boxplot()+
  theme_classic2()+
  scale_color_manual(values = c("#65a86c", "#a865a1" ))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=8),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 9))+
  labs( x = "Genotype", y = "log tlc1 to CaM ratio",
       color="Ratio Change")


#combined plot

full_facet <- ggarrange(facet_g, full_g, font.label = list(size = 10, color = "black"), labels=c("A", "B", "C"), common.legend = FALSE)

ggsave("/Volumes/RSMASFILES2/ThesisPlots/fig_2.png", width = 11, height = 8, units = "in",
       dpi=300)



###
qPCR_16S <- lm(rica_16S_abund~TLC_CAM_log, data= rica)
summary(qPCR_16S)

rica1 <- rica %>% filter(TLC_CAM_ratio<2)
ggplot(rica, aes(x=rica_16S_abund, y=TLC_CAM_log ))+
  geom_point()+
  stat_smooth(method = lm, color="cyan")+
  stat_regline_equation(aes(label =  paste( ..adj.rr.label.., sep = "~~~~")))+
  theme_classic2()+
  labs(y="log tlc1 to CaM ratio", x="16S MD3-55 Abundance", 
       caption = "Multiple R-squared:  0.217,  Adjusted R-squared:  0.2081 
F-statistic: 24.39 on 1 and 88 DF,  p-value: 3.711e-06")+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 12),
        plot.caption=element_text(size = 10))

ggsave("/Volumes/RSMASFILES2/ThesisPlots/16s_qPCRmodel.png", width = 5.5, height = 5, units = "in",
       dpi=300)





###
##dPCR 
###

rica1 <- rica %>% filter(Time_Point=="T2") 

dPCR_16S <- lm(rica_16S_abund~dPCR_conc_thermo, data= rica1)
summary(dPCR_16S)

ggplot(rica1, aes(x=rica_16S_abund, y=dPCR_conc_thermo ))+
  geom_point()+
  stat_smooth(method = lm, color="cyan")+
  stat_regline_equation(aes(label =  paste( ..adj.rr.label.., sep = "~~~~")))+
  theme_classic2()+
  labs(y="dPCR tlc1 absolute adundance", x="16S MD3-55 Abundance", 
       caption = "Multiple R-squared:  0.05212, Adjusted R-squared:  0.0265 
F-statistic: 2.034 on 1 and 37 DF,  p-value: 0.1622")+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 12),
        plot.caption=element_text(size = 10))

