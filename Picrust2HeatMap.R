library(ggpicrust2)
library(tidyverse)
library(readr)
library(ggprism)
library(patchwork)
library(ggh4x)
library(phyloseq)
library(microViz)
library(ComplexHeatmap)

####
## Picrust2 analysis for RICA data (qiime2 processed)
## NOAA - AOML - University of Miami
## Sterling Butler
## 06042024
####

#look at the data supplied by ggpicrust2
cyc<-data(metacyc_abundance)
ko <-data(ko_abundance)
data(metadata)


# Set directory 
setwd("/Users/sterling.butler/Desktop/RICA_picrust2_files/picrust2_out_pipeline/pathways_out")

#Get abundance file
abundance_file <- read_delim("path_abun_unstrat.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE) 

#Get metadata
meta_df<- read.csv("qiime203112024_16S_RICA_META.csv")

# Set directory 
setwd("/Users/sterling.butler/Desktop/RICA_picrust2_files/picrust2_out_pipeline/EC_metagenome_out")
abundance_file_pred <- read_delim("pred_metagenome_unstrat_descrip.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)

#get reference pathway map file
load("/Users/sterling.butler/Downloads/MetaCyc_pathway_map.RData")
pathway_map <- MetaCyc_pathway_map

#make the correct column names for the phyloseq object 
abundance_file1 = abundance_file %>% tibble::column_to_rownames("pathway")
meta = meta_df %>% tibble::column_to_rownames("sample")

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

#make a pcOA with time point
phy_clr = microbiome::transform(pi, 'clr')
phy_ord= ordinate(phy_clr, "RDA", "euclidean")
plot_ordination(pi, phy_ord, color= "Spirochaetaceae_ANCOMpres")

#make a NMDS with time point
phy_ord2= ordinate(phy_clr, "NMDS", "euclidean")
plot_ordination(pi, phy_ord2, color= "Time_Point", )

#plot_ordination##test if they are different
library("pairwiseAdonis")
dist.uf <- phyloseq::distance(phy_clr, method = "euclidean")
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Time_Point, sim.method = "euclidean",
                p.adjust.m = "bonferroni")


?pairwise.adonis2
#heat map for pathways
pi_phylum = tax_glom(pi, taxrank = "Superclass2",NArm=FALSE)
plot_heatmap(pi_phylum,"NMDS", "bray", "Time_Point", "Superclass2", first.sample="T1")+
  facet_grid(.~Time_Point, scales = "free")


########
#ANCOM-BC2 with just timepoint no grouping and using ssuperclass not by pathways  
########
library(ANCOMBC)
pi_ancom = mia::makeTreeSummarizedExperimentFromPhyloseq(pi)
print(pi_ancom)


set.seed(123)

output3 = ancombc2(
  data = pi_ancom,
  assay_name = "counts",
  tax_level = "Superclass2",
  fix_formula = "Time_Point",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Time_Point",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = FALSE,  # Deactivated due to <3 categories
  pairwise = FALSE,  # Deactivated due to <3 categories
  dunnet = FALSE,  # Deactivated due to <3 categories
  trend = FALSE,  # Deactivated due to <3 categories
  iter_control = list(
    tol = 1e-2,
    max_iter = 20,
    verbose = TRUE
  ),
  em_control = list(
    tol = 1e-5,
    max_iter = 100
  ),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(
    fwer_ctrl_method = "holm",
    B = 100))

res_prim3 = output3$res

df_fig_TP3 = res_prim3 %>%
  dplyr::filter(diff_Time_PointT2 == "TRUE") %>% 
  dplyr::arrange(desc(lfc_Time_PointT2)) %>%
  dplyr::mutate(direct = ifelse(lfc_Time_PointT2 > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_Time_PointT2 == "TRUE", "aquamarine3", "black"))


colnames(pathway_map)[1]<- "taxon"
df_fig_TP3 <- left_join(df_fig_TP3, pathway_map, by="taxon")
df_fig_TP3$taxon = factor(df_fig_TP3$taxon, levels = df_fig_TP3$taxon)


lfc_path <-df_fig_TP3 %>% ggplot(aes(x = taxon, y = lfc_Time_PointT2, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin = lfc_Time_PointT2 - se_Time_PointT2, ymax = lfc_Time_PointT2 + se_Time_PointT2), 
                width = 0.2, position = position_dodge(0.05), color = "black")+
  coord_flip()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 15),
        title = element_text(size = 15))+
  labs(y="Log Fold Change", x="Pathway", title = "Log Fold Change from T1 to T2")



#Volcano plot at the pathway level

#create an ancomc for pathway level to make volcano 
set.seed(123)

output2 = ancombc2(
  data = pi_ancom,
  assay_name = "counts",
  tax_level = "pathway",
  fix_formula = "Time_Point",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Time_Point",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = FALSE,  # Deactivated due to <3 categories
  pairwise = FALSE,  # Deactivated due to <3 categories
  dunnet = FALSE,  # Deactivated due to <3 categories
  trend = FALSE,  # Deactivated due to <3 categories
  iter_control = list(
    tol = 1e-2,
    max_iter = 20,
    verbose = TRUE
  ),
  em_control = list(
    tol = 1e-5,
    max_iter = 100
  ),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(
    fwer_ctrl_method = "holm",
    B = 100))

res_prim2 = output2$res


#label the data 
res_prim2_lab <- res_prim2 %>%
  mutate(reg_status = case_when(lfc_Time_PointT2 < (-0.5)  ~ "Down Regulated",
                                lfc_Time_PointT2 > 0.5  ~ "Up Regulated",))


res_prim3_lab$delabel <- ifelse(res_prim3_lab$taxon %in% head(res_prim3_lab[order(res_prim3_lab$p_Time_PointT2), "taxon"], 9), res_prim3_lab$taxon, NA)

#create the plot 
path_volc <- ggplot(res_prim2_lab, aes(x=lfc_Time_PointT2, y= -log10(p_Time_PointT2), col = reg_status))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed')+ 
  scale_color_manual(values = c("#00AFBB", "#bb0c00", "grey"), # to set the colours of our variable  
                     labels = c("Decrease", "Increase" , "Not significant"))+
  labs(title = "Volcano Plot for Picrust2 data", 
       x = "Log Fold Change from T1 to T2", y = "-log10(p-value)",
       color="Pathway Abundance Status")+
  theme_classic2()+
  theme(axis.text.x = element_text( size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 15),
        title = element_text(size = 15))


#combined the two graphs 

ggarrange(lfc_path , path_volc, font.label = list(size = 10, color = "black"), labels=c("A", "B"), common.legend = FALSE)



#####
## same same but for the superclass level used to show lfc graph
####
#label the data 
res_prim3_lab <- res_prim3 %>%
  mutate(reg_status = case_when(lfc_Time_PointT2 < (-0.5) ~ "Down Regulated",
                                lfc_Time_PointT2 > 0.5  ~ "Up Regulated",))


res_prim3_lab$delabel <- ifelse(res_prim3_lab$taxon %in% head(res_prim3_lab[order(res_prim3_lab$p_Time_PointT2), "taxon"], 9), res_prim3_lab$taxon, NA)

#create the plot 
ggplot(res_prim3_lab, aes(x=lfc_Time_PointT2, y= -log10(p_Time_PointT2), col = reg_status))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed')+ 
  scale_color_manual(values = c("#00AFBB", "#bb0c00", "grey"), # to set the colours of our variable  
                     labels = c("Downregulated", "Upregulated" , "Not significant"))+
  labs(title = "Volcano Plot for Picrust2 data", 
       x = "Log Fold Change from T1 to T2", y = "-log10(p-value)",
       color="Pathway Status")



#pathway abundance 


library("randomForestSRC")

predictors <- t(otu_table(pi))
response <- as.factor(sample_data(pi)$Time_Point)
rf.data <- data.frame(response, predictors)
head(rf.data, n=2)


tp.classify <- rfsrc(response~., data = rf.data, ntree = 2000,
                       importance="permute", csv.num=TRUE)


tp.classify
plot(tp.classify)

ggsave("/Volumes/RSMASFILES2/ThesisPlots/tp.classify.png", width = 11, height = 9, units = "in",
       dpi=300)

##Combine top important ASVs from each reef
set.seed(2)
# Create a list to store the results
result_list <- list()

# Define the column names for sorting
sort_columns <- c("T1", "T2")

# Loop through the sorting columns
for (col in sort_columns) {
  result <- cbind(as.data.frame(tp.classify$importance), as.data.frame(tax_table(pi))) %>%
    arrange(desc(get(col))) %>%
    head(n=5) %>%
    distinct()  # Remove duplicate rows
  
  # Assign the result to the list
  result_list[[paste0(col, "_5")]] <- result
}

# Combine the results into a single data frame
tp_df <- do.call(rbind, result_list)

rf_trip_df <- tp_df %>%
  rownames_to_column("ID") %>%
  separate("ID", into = c("Time_Point", "pathway"), sep = "\\.", remove = FALSE) %>% 
  unite("Trip_top", Time_Point, na.rm = TRUE, remove = FALSE) %>%
  column_to_rownames("ID")
 
df_trip_long<-result  %>%
  rownames_to_column() %>%
  pivot_longer(cols = c(T1, T2), names_to = "variable") %>%
  group_by(variable)


#make a plot of importance valued all overall imporatance 

ggplot(df_trip_long,
  aes(x=all, y = value)) + 
  geom_point(aes(shape=variable, color=pathway ), 
             alpha=0.7,
             size=4) +
  scale_color_manual("Pathway", 
      values=c( "#AD6F3B", "#6F8FAF", "#800020", "lightpink", "#56B4E9")) +
  theme_classic2()+
  ylab("Pathway Importance per Time Point") +
  xlab("Overall Pathway Importance") 

ggsave("/Volumes/RSMASFILES2/ThesisPlots/rf_picrust.png", width = 11, height = 9, units = "in",
       dpi=300)



