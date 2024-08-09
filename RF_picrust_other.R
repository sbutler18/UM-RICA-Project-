library(ggpicrust2)
library(tidyverse)
library(readr)
library(ggprism)
library(patchwork)
library(ggh4x)
library(phyloseq)




# Set directory 
setwd("/Users/sterling.butler/Desktop/RICA_picrust2_files/picrust2_out_pipeline/pathways_out")

#Get metadata
meta_df<- read.csv("qiime203112024_16S_RICA_META.csv")

# Set directory 
setwd("/Users/sterling.butler/Desktop/RICA_picrust2_files/picrust2_out_pipeline/EC_metagenome_out")
abundance_file <- read_delim("pred_metagenome_unstrat_descrip.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_file$`function` <- NULL


#get reference pathway map file
load("/Users/sterling.butler/Downloads/MetaCyc_pathway_map.RData")
row.names(MetaCyc_pathway_map) <- NULL

#make the correct column names for the phyloseq object 
abundance_file1 <- abundance_file %>%
  mutate(description = make.unique(description)) %>%
  tibble::column_to_rownames("description") 
  

meta = meta_df %>% tibble::column_to_rownames("sample")

pathway_map <- MetaCyc_pathway_map %>% 
  mutate(description = make.unique(pathway)) %>%
  tibble::column_to_rownames("pathway")


#convert the taxa table back into a matrix, so phyloseq can read
abun_mat = as.matrix(abundance_file1)
pathway_map_mat = as.matrix(pathway_map)

#make the individual phyloseq objects otu_table, tax_table, sample_data
otu<-otu_table(abun_mat, taxa_are_rows = TRUE)
tax_table<-tax_table(pathway_map_mat)
sample<-sample_data(meta)

#Check to see if the labels for samples and ASVs for the different 
sample_names(sample)
sample_names(otu)
taxa_names(tax_table)
taxa_names(otu)

#lets make that phyloseq object now
pi= phyloseq(otu, tax_table, sample)

#filter
pi <-pi %>% subset_samples(Species!="NA")
pi <-pi %>% subset_samples(Sample.NTC.Control!="32")



#pathway abundance 


library("randomForestSRC")

predictors <-as.data.frame(t(abundance_file1))
response <- as.factor(meta$Time_Point)
rf.data <- data.frame(response, predictors)
head(rf.data, n=2)


tp.classify <- rfsrc(response~., data = rf.data, ntree = 3000,
                     importance="permute", csv.num=TRUE)


tp.classify
plot(tp.classify)
