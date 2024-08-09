library(phyloseq)
library(ggplot2)
library(qiime2R)
library(microViz)
library(tidyverse)
library(ggpubr)


###
# Sterling Butler
# NOAA/ UM 
# RICA Nursery Project - qiime2 processed barplots with microViz
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
 select(-Confidence) %>% 
 arrange(Feature.ID) %>% 
 mutate(ASV = 1:n()) %>% 
 mutate(newcol = "ASV") %>%
 unite("ASVs", newcol:ASV)

#make the taxa table into a data frame 
taxa<-as.data.frame(taxa)


setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged")

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

#remove NTC and a sample with 72 total reads 
phy <-phy %>% subset_samples(Species!="NA")
phy <-phy %>% subset_samples(Bag_Number!="32")

#filter out the noise in the data, this filters out any taxa that have counts less than 5
phy_fill= filter_taxa(phy , function(x) sum(x > 5) > (0.05*length(x)), TRUE)

## Basic summary information for the phyloseq object that you just created
#number of taxa
ntaxa(phy_fill)
ntaxa(phy)

#number of samples
nsamples(phy_fill)
#number of variables
sample_variables(phy_fill)

#Use tax_fix() on your phyloseq data with default arguments to repair most tax_table problems (missing or uninformative values)
phy_fill <- tax_fix(phy_fill)

phy_fill 


#Relative Abundance for T1 and T2 
phy_fill %>%
  comp_barplot(tax_level = "Family",
               label = "Genotype",
               other_name = "Other genera") +
  facet_grid(~Time_Point, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))+
  labs(x = "Genotype", y = "Relative Abundance")+
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent())



#T1 and T2 facet by genotype
phy_fill %>% 
  ps_arrange(Genotype, Time_Point)%>%
  ps_mutate(Time_Point=as.numeric(Time_Point)) %>%
  comp_barplot(tax_level = "Family",
               label = "Time_Point",
               other_name = "Other genera"
               ) +
  facet_grid(~Genotype, scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))+
  labs(x = "Time_Point", y = "Relative Abundance")+
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent())



#Not using the microViz comp_barplot
df_test   <- phy_fill %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()

head(df_test)

#all taxa for family

all_tax<- df_test %>%
  filter(Abundance >0.005) %>%
  ggplot(aes(x =Time_Point, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity", position="fill", aes(fill = Family)) + 
  facet_wrap(.~a_Genotype, scales = "free") +
  scale_fill_manual(values=c( "#56B4E9","#CBD588", "orange","#DA5724","#CD9BCD",
    "gray80", "#AD6F3B", "green","#D14285", "#652926","#8569D5", 
    "#5E738F","#D1A33D", "#8A7C64","lightsalmon","aquamarine4",
    "lightblue4", "lightpink", "ivory4","royalblue4", "darkorchid", 
    "palevioletred1", "#56B4E9","#CBD588","yellow2","#5F7FC7", "orange","#DA5724",
    "#CD9BCD", "gray80",
    "#AD6F3B" ),
  labels = c("Unknown Bacteria", "Fastidiosibacter lacustris", "Aquarickettsia (MD3-55)",
             "Candidatus Gracilibacteria", "Rhodospirillaceae", "Roseobacteraceae", 
             "SAR324 clade (Marine_group_B)", "Spirochaetaceae", "Thalasssobaculaceae",
             "Thiovulaceae", "Vibronaceae", "Campylobacterales", "Unknown Bacteria")) +
  guides(fill = guide_legend(keywidth = 0.3, , keyheight =.50, ncol=1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
        strip.text.x = element_text(
          size = 9, face = "bold.italic")) +
  theme(strip.background = element_rect(
    fill="white", color="black", size=.5))+
  theme(legend.text = element_text(size=12)) +
  labs(color = "Bacteria", title = "Relative Adundance (All Species)", x= "Time Point")



ggsave("/Volumes/RSMASFILES2/ThesisPlots/fig_4.png", width = 9, height = 8, units = "in",
       dpi=300)


#No Aquarickettsia
rica_no_tax <- df_test %>%
filter(Family!=" f__Fokiniaceae") %>%
filter(Abundance >0.005) %>%
ggplot(aes(x =Time_Point, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity", position="fill", aes(fill = Family)) + 
  facet_wrap(.~a_Genotype, scales = "free") +
  scale_fill_manual(values=c("#56B4E9","#CBD588","#DA5724","#CD9BCD",
                             "gray80", "#AD6F3B", "green","#D14285", "#652926","#8569D5", 
                             "#5E738F","#D1A33D", "#8A7C64","lightsalmon","aquamarine4",
                             "lightblue4", "lightpink", "ivory4","royalblue4", "darkorchid", 
                             "palevioletred1", "#56B4E9","#CBD588","yellow2","#5F7FC7", "orange","#DA5724",
                             "#CD9BCD", "gray80",
                             "#AD6F3B"),
                            labels = c("Unknown Bacteria", "Fastidiosibacter lacustris",
                                      "Candidatus Gracilibacteria", "Rhodospirillaceae", "Roseobacteraceae", 
                                     "SAR324 clade (Marine_group_B)", "Spirochaetaceae", "Thalasssobaculaceae",
                                       "Thiovulaceae", "Vibronaceae", "Campylobacterales", "Unknown Bacteria"
                             )) +
  guides(fill = guide_legend(keywidth = 0.3, , keyheight =.50, ncol=1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
        strip.text.x = element_text(
          size = 10, face = "bold.italic")) +
  theme(strip.background = element_rect(
    fill="white", color="black", size=.5))+
  theme(legend.text = element_text(size=12)) +
  labs(color = "Bacteria", title = "Relative Adundance (No Aquarickettsia)", x= "Time Point")


ggsave("/Volumes/RSMASFILES2/ThesisPlots/fig_5.png", width = 9, height = 8, units = "in",
       dpi=300)


#combine both graphs

a_b <- ggarrange(all_tax, rica_no_tax, font.label = list(size = 10, color = "black"), labels=c("A", "B"), common.legend = FALSE)


ggsave("/Volumes/RSMASFILES2/ThesisPlots/fig_4.png", width = 8, height = 4.5, units = "in",
       dpi=300)

