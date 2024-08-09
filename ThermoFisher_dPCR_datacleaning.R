
###
#NOAA - Sterling Butler
#dPCR pipeline
#takes thermofisher dPCR outputs and combines multiple plates into one document 
#08/2024
###


library(tidyverse)
library("dplyr")                                                 
library("plyr")                                                
library("readr")


#upload and combined all plates
#set path to location of thermofisher files
dPCR_results <- list.files(path="/Users/sterling.butler/Desktop/dPCR_Data/thermofisher", pattern=".csv", full.names=T)%>%
  lapply(read_csv) %>%
  bind_rows()

#set column names
dPCR_results <- dPCR_results %>% set_names(dPCR_results[1, ])

#remove other column names from dataset from different plates
dPCR_results <- dPCR_results %>% filter(Well!="Well")

names(dPCR_results)
#check NTCs

# 1. Check and remove NTC wells
ntc <- dPCR_results[which(dPCR_results$Sample=="NTC"), ]

dPCR_results <- dPCR_results %>% filter(Sample!="NTC")
dPCR_results <- dPCR_results %>% filter(Sample!="ntc")
# 2. Investigate blank values within the dataframe
blank <- dPCR_results[!complete.cases(dPCR_results$Conc.),  ]

dPCR_results <- dPCR_results %>% filter(Conc.!="NA")

# 3. Investigate values that were out side of a set Precision % level, here is set to select any run that have >5% percision 
thres <- 5
sketch <- dPCR_results %>% dplyr::filter( `Precision %`  < thres)

dPCR_results1 <- dPCR_results %>%  dplyr::filter(`Precision %` > 2)



write.csv(dPCR_results, "/Users/sterling.butler/Desktop/dPCR_Data/thermofisher/allplate_rica.csv",
    )





