library(tidyverse) ; packageVersion("tidyverse") # 1.3.2
library(phyloseq) ; packageVersion("phyloseq") # 1.42.0
library(ggpubr)
library(microbiome)
library(ggplot2)
library(qiime2R)
library(microViz)

devtools::install_github("AndreaCirilloAC/updateR")
updateR()

###
# Sterling Butler
# NOAA/ UM 
# RICA Nursery Project - qiime2 processed PCOA plots with microViz
# 3/2024
###

#set working directory
setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged/feature")

#import the feature table using the read_q2biom()
asv<- read_q2biom("table.biom")

#import the taxa file as a qiime2 artifact and then convert it into a usable tax table for phyloseq
setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged")
tax <- read_qza("taxonomy.qza")
taxa <- tax$data %>% as_tibble() %>% 
  separate(Taxon, sep=";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  select(-Confidence) %>% arrange(Feature.ID) %>% mutate(ASV = 1:n()) %>% 
  mutate(newcol = "ASV") %>%
  unite("ASVs", newcol:ASV)

#make the taxa table into a data frame 
taxa<-as.data.frame(taxa)

#import the meta data
meta_df<-read.csv("qiime2_08022024_16S_RICA_META.csv")

#make the correct column names for the phyloseq object 
taxa1 = taxa %>% tibble::column_to_rownames("Feature.ID")
meta = meta_df %>% tibble::column_to_rownames("sample")

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

#Use tax_fix() on your phyloseq data with default arguments to repair most tax_table problems (missing or uninformative values)
phy_fill <- tax_fix(phy_fill)

#####standard PCA for t1 
phy_fill %>% 
  subset_samples(Time_Point=="T1") %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  ord_calc("PCA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 10))

#ordiations plot with ord_calc for standard PCA
phy_fill %>%
  subset_samples(Time_Point=="T1") %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2),color = "Genotype", plot_taxa = 1:3, size = 2.5,
           tax_lab_style = tax_lab_style(size = 3, alpha = 0.5))+
  scale_colour_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "darkorchid","#DA5724","#CD9BCD","#006700","yellow2",
                               "gray80", "#AD6F3B", "#673770","#8569D5", 
                               "#5E738F","#D1A33D", "orange","#ff6700","aquamarine4", "#652926",
                               "lightblue4", "lightpink","royalblue4","#D14285",
                               "palevioletred1", "#56B4E9","#CBD588", "#5F7FC7","#DA5724",
                               "#CD9BCD", "gray80", "darkorchid",
                               "#AD6F3B", "#673770","#D14285", "#652926"))


#####standard PCA for t1 to t2
#look at distibution of the eigenvalues, for this data most was in PC1-PC2, so I am selecting those by using axes = c(1, 2) when plotting below
phy_fill %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  ord_calc("PCA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 10))

#ordiations plot with ord_calc for standard PCA
phy_fill %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2),color = "Spirochaetaceae_ANCOMpres", shape="Time_Point",  size = 2.5,
           tax_lab_style = tax_lab_style(size = 4, alpha = 0.5))+
  scale_colour_manual(values=c("#000000", "#7de661"))+
  theme(axis.text.x = element_text( size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 15),
        title = element_text(size = 15))+
  labs(title = "Principal Component Analysis for Timepoints and Spirochaetaceae", color= "Spirochaetaceae enriched", shape= "Time Point")
  
ggsave("/Volumes/RSMASFILES2/ThesisPlots/spiro_pcoa.png", width = 11, height = 7, units = "in",
       dpi=300)  

#test the if groups are different 

#Time Point 
## pvalue=0.001
phy_clr = microbiome::transform(phy_fill, 'clr')
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Time_Point, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#Spirochaetaceae_ANCOMpres
## pvalue=0.001
phy_clr = microbiome::transform(phy_fill, 'clr')
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Spirochaetaceae_ANCOMpres, sim.method = "euclidean",
                p.adjust.m = "bonferroni")


phy_fill %>%
  tax_transform("clr", rank = "Genus") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2),color = "Genotype", shape="Time_Point", plot_taxa = 1:3, size = 2.9,
           tax_lab_style = tax_lab_style(size = 4, alpha = 0.5))+
  scale_colour_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "darkorchid","#DA5724","#CD9BCD","#006700","yellow2",
                               "gray80", "#AD6F3B", "#673770","#8569D5", 
                               "#5E738F","#D1A33D", "orange","#ff6700","aquamarine4", "#652926",
                               "lightblue4", "lightpink","royalblue4","#D14285",
                               "palevioletred1", "#56B4E9","#CBD588", "#5F7FC7","#DA5724",
                               "#CD9BCD", "gray80", "darkorchid",
                               "#AD6F3B", "#673770","#D14285", "#652926"))+
  theme_classic2()+
  theme(axis.text.x = element_text( size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 15),
        title = element_text(size = 15))+
  labs(title = "Principal Component Analysis for Genotypes and Timepoints", shape= "Time Point")

######Aitchison Distance 
#look at distibution of the eigenvalues
phy_fill %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6))

#ordiations plot with ord_calc, Aitchison Distance 
phy_fill %>%
    tax_transform("identity", rank = "Genus") %>% # don't transform!
    dist_calc("aitchison") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(color  = "Time_Point", plot_taxa = 1:5, size = 2.9)



