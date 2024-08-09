library(tidyverse)


#set working directory
setwd("/Volumes/RSMASFILES2/RICA_DATA_NEW/")

#get metadata
rica <- read.csv("06262024_RICA_Nursery_qPCR_META.csv")
head(rica)

#high means 
rica_high_16S <- rica %>% filter(Genotype=="Acerv-1" | Genotype== "Marker 9" | Genotype== "Steph's A - B4"
                                 |Genotype=="Acerv-3"|Genotype=="Cooper's"| Genotype=="Sunny Isles F"|
                                   Genotype=="Cheetos D")

mean_high <- round(mean(rica_high_16S$TLC_CAM_ratio), 3)
sd_high <- round(sd(rica_high_16S$TLC_CAM_ratio),3)

#low means
rica_low_16S <- rica %>% filter(Genotype== "Sunny Isles A"|Genotype== "Acerv-2" |Genotype== "Thicket-3"|
                                  Genotype== "Sunny Isles C"|Genotype== "Sunny Isles D"|
                                  Genotype== "Sandy" |Genotype== "Tuna Jelly-2"|Genotype== "Stag Reef A")

mean_low <- round(mean(rica_low_16S$TLC_CAM_ratio), 3)
sd_low <- round(sd(rica_low_16S$TLC_CAM_ratio), 3)

#mid means 
rica_mid_16S <- rica %>% filter(Genotype== "KJS-74" | Genotype== "Elkhorn" | Genotype== "Miami Beach B" | 
                                  Genotype== "Steph's A - C6" | Genotype== "Stag C"| 
                                  Genotype== "Kaufman-1" | Genotype== "Kelsey's 1" |Genotype== "Cheetos A" )

mean_mid <- round(mean(rica_mid_16S$TLC_CAM_ratio), 3)
sd_mid <- round(sd(rica_mid_16S$TLC_CAM_ratio), 3)