library(tidyverse) ; packageVersion("tidyverse") # 1.3.2
library(phyloseq) ; packageVersion("phyloseq") # 1.42.0
library(ggpubr)
library(microbiome)
library(ggplot2)
library(qiime2R)

#######
## Sterling Butler 
## NOAA, RSMAS University of Miami 
## 12/14/2023
## Aquarickettsia x Acropora cervicornis  
## UM nursery models and statistical anaylsis
#######

#set working directory
setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged/feature")

#import the feature table using the read_q2biom()
asv<- read_q2biom("table.biom")

#import the taxa file as a qiime2 artifact and then convert it into a usable tax table for phyloseq
setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged")
tax <- read_qza("taxonomy.qza")
taxa <- tax$data %>% as_tibble() %>% 
  separate(Taxon, sep=";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  select(-Confidence) %>% 
  arrange(Feature.ID) %>% 
  mutate(ASV = 1:n()) %>% 
  mutate(newcol = "ASV") %>%
  unite("ASVs", newcol:ASV)

#make the taxa table into a data frame 
taxa<-as.data.frame(taxa)

#import the meta data
meta_df<-read.csv("qiime2_05202024_16S_RICA_META.csv")


#make the correct column names for the phyloseq object 
taxa1 = taxa %>% tibble::column_to_rownames("Feature.ID")
meta = meta_df %>% tibble::column_to_rownames("sample")

#factor the t1/t2
meta$Time_Point <- factor(meta$Time_Point, levels = c("T1","T2"))

#convert the taxa table back into a matrix, so phyloseq can read
tax_mat = as.matrix(taxa1)

#make the individual phyloseq objects otu_table, tax_table, sample_data
otu<-otu_table(asv, taxa_are_rows = TRUE)
tax_table<-tax_table(tax_mat)
sample<-sample_data(meta)

#Check to see if the labels for samples and ASVs for the different 
sample_names(sample)
sample_names(otu)
taxa_names(tax_table)
taxa_names(otu)

#lets make that phyloseq object now
phy= phyloseq(otu, tax_table, sample)

#remove Mitochondria/ chloroplast
phy <- phy %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 

#remove NTC and a sample with 20 total reads 
phy <-phy %>% subset_samples(Species!="NA")
phy <-phy %>% subset_samples(Bag_Number!="32")

#filter out the noise in the data, this filters out any taxa that have counts less than 5
phy_fill= filter_taxa(phy , function(x) sum(x > 5) > (0.05*length(x)), TRUE)

## Basic summary information for the phyloseq object that you just created
#number of taxa
ntaxa(phy_fill)
#number of samples
nsamples(phy_fill)
#number of variables
sample_variables(phy_fill)


#Spirochateaece 
#abundace for aquarickettsia spiro yes
Spiro_yes <- subset_samples(phy_fill, Spirochaetaceae_ANCOMpres == "NO")
Spiro_yes_psmelt <- Spiro_yes %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#yes mean/sd
rica_spiro_yes <- Spiro_yes_psmelt %>% filter(OTU=="bd53e2483c985a090089d93fc039c4ac") 
mean_rica_spiro_yes <-  round(mean(rica_spiro_yes$Abundance), 3)
sd_rica_spiro_yes <-  round(sd(rica_spiro_yes$Abundance), 3)

#abundace for aquarickettsia spiro no
Spiro_no <- subset_samples(phy_fill, Spirochaetaceae_ANCOMpres == "YES")
Spiro_no_psmelt <- Spiro_no %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#no mean/sd
rica_spiro_no <- Spiro_no_psmelt %>% filter(OTU=="bd53e2483c985a090089d93fc039c4ac") 
mean_rica_spiro_no <-  round(mean(rica_spiro_no$Abundance), 3)
sd_rica_spiro_no <-  round(sd(rica_spiro_no$Abundance), 3)


#T1 to T2 relative abundance

#t1
t1_ra <- subset_samples(phy_fill, Time_Point == "T1")
t1_ra_psmelt <- t1_ra %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#t1 campylobacterales mean/sd
campylo_t1_ra <- t1_ra_psmelt %>% filter(OTU=="70d69f68d0230abfbb2f096379d6131e") 
mean_campylo_t1_ra <-  round(mean(campylo_t1_ra$Abundance), 3)
sd_campylo_t1_ra <-  round(sd(campylo_t1_ra$Abundance), 3)

#t2
t2_ra <- subset_samples(phy_fill, Time_Point == "T2")
t2_ra_psmelt <- t2_ra %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#t2 campylobacterales mean/sd
campylo_t2_ra <- t2_ra_psmelt %>% filter(OTU=="70d69f68d0230abfbb2f096379d6131e") 
mean_campylo_t2_ra <-  round(mean(campylo_t2_ra$Abundance), 3)
sd_campylo_t2_ra <-  round(sd(campylo_t2_ra$Abundance), 3)


###
#Rhodospirillaceae t1 and t2
###

view(tax_table(phy_fill))


rhodo_t1 <-  t1_ra %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#t1 mean/sd
rhodo_t1 <- rhodo_t1 %>% filter(OTU=="5bbc92a4955ba130b3d4157df25a2608") 
mean_rhodo <-  round(mean(rhodo_t1$Abundance), 3)
sd_rhodo <-  round(sd(rhodo_t1$Abundance), 3)


rhodo_t2 <-  t2_ra %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#t2 mean/sd
rhodo_t2 <- rhodo_t2 %>% filter(OTU=="5bbc92a4955ba130b3d4157df25a2608") 
mean_rhodo_2 <-  round(mean(rhodo_t2$Abundance), 3)
sd_rhodo_2 <-  round(sd(rhodo_t2$Abundance), 3)



#Spirochateaece
#t1 mean/sd
spiro_t1 <-  t1_ra %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
spiro_t1 <- spiro_t1 %>% filter(OTU=="3aa965194c47dee2ebb9ad04c942930a") 
mean_spiro_1 <-  round(mean(spiro_t1$Abundance), 3)
sd_spiro_1 <-  round(sd(spiro_t1$Abundance), 3)

#t2 mean/sd
spiro_t2 <-  t2_ra %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
spiro_t2 <- spiro_t2 %>% filter(OTU=="3aa965194c47dee2ebb9ad04c942930a") 
mean_spiro_2 <-  round(mean(spiro_t2$Abundance), 3)
sd_spiro_2 <-  round(sd(spiro_t2$Abundance), 3)



#RICA Threshold abundance 
#set at 70% 
rica_abund <- phy_fill %>%  tax_glom(taxrank = "Order") %>%  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

rica_abund2 <- rica_abund %>% filter(OTU == c( 
  "bd53e2483c985a090089d93fc039c4ac"))

rica_abund_filt <- rica_abund2 %>% select(Abundance, Bag_Number)
rica_new<-left_join(meta_df, rica_abund_filt, by="Bag_Number")

rica_new <- rica_new %>%
  mutate(rica_ab_level = if_else(Abundance>0.70, "HIGH", "LOW"))

meta2 = rica_new %>% tibble::column_to_rownames("sample")
sample2<-sample_data(meta2)


#lets make that phyloseq object now
phy2= phyloseq(otu, tax_table, sample2)

#remove Mitochondria/ chloroplast
phy2 <- phy2 %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 

#remove NTC and a sample with 20 total reads 
phy2 <-phy2 %>% subset_samples(Species!="NA")
phy2 <-phy2 %>% subset_samples(Bag_Number!="32")

#filter out the noise in the data, this filters out any taxa that have counts less than 5
phy_fill2= filter_taxa(phy2 , function(x) sum(x > 5) > (0.05*length(x)), TRUE)

#low rica adundance tax 
low_rica <- subset_samples(phy_fill2, rica_ab_level == "LOW")
low_rica_psmelt <- low_rica %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#thassobacter low rica mean/sd
Thasso_low <- low_rica_psmelt %>% filter(OTU=="d4a521b98cbc2f7183b70400d38c0216") 
mean_Thasso_low <-  round(mean(Thasso_low$Abundance), 3)
sd_Thasso_low <-  round(sd(Thasso_low$Abundance), 3)

#high rica adundance tax 
high_rica <- subset_samples(phy_fill2, rica_ab_level == "HIGH")
high_rica_psmelt <- high_rica %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#no mean/sd
thasso_high <- high_rica_psmelt %>% filter(OTU=="d4a521b98cbc2f7183b70400d38c0216") 
mean_thasso_high <-  round(mean(thasso_high$Abundance), 3)
sd_thasso_high <-  round(sd(thasso_high$Abundance), 3)


###
#Rhodospirillaceae
###

view(tax_table(phy_fill))


rhodo_t1 <-  t1_ra %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

#yes mean/sd
rhodo_t1 <- rhodo_t1 %>% filter(OTU=="5bbc92a4955ba130b3d4157df25a2608") 
mean_rhodo <-  round(mean(rhodo_t1$Abundance), 3)
sd_rhodo <-  round(sd(rhodo_t1$Abundance), 3)


rhodo_t2 <-  t2_ra %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#yes mean/sd
rhodo_t2 <- rhodo_t2 %>% filter(OTU=="5bbc92a4955ba130b3d4157df25a2608") 
mean_rhodo_2 <-  round(mean(rhodo_t2$Abundance), 3)
sd_rhodo_2 <-  round(sd(rhodo_t2$Abundance), 3)


