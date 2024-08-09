library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(DT)

###
# Sterling Butler
# NOAA/ UM 
# RICA Nursery Project - qiime2 processed ANCOM-BC2 
# 5/2024
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
meta_df<-read.csv("qiime2_05202024_16S_RICA_META.csv")

#make the correct column names for the phyloseq object 
taxa1 = taxa %>% tibble::column_to_rownames("Feature.ID")
meta = meta_df %>% tibble::column_to_rownames("sample")

#convert the taxa table back into a matrix, so phyloseq can read
tax_mat = as.matrix(taxa1)

#make the individual phyloseq objects otu_table, tax_table, sample_data
otu<-otu_table(asv, taxa_are_rows = TRUE)
tax_table<-tax_table(tax_mat)
sample<-sample_data(meta)

#lets make that phyloseq object now
phy= phyloseq(otu, tax_table, sample)

#remove Mitochondria/ chloroplast
phy <- phy %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 
## Basic summary information for the phyloseq object that you just created
#number of taxa
ntaxa(phy)
#number of samples
nsamples(phy)
#number of variables
sample_variables(phy)



#remove NTC and a sample with 20 total reads 
#phy <-phy %>% subset_samples(Species!="NA")
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

#make into treesummarized fromat for ancom function
rica_ancom = mia::makeTreeSummarizedExperimentFromPhyloseq(phy)
print(rica_ancom)

# Subset to genotypes
rica_ancom_g = rica_ancom[, rica_ancom$a_Genotype %in% c("Acerv-1", "Acerv-2", "Acerv-3", "Cheetos A", "Cheetos D", "Cooper's",
                                                       "Elkhorn", "Kaufman-1", "Kelsey's 1", "KJS-74", "Marker 9",
                                                       "MB B", "Sandy", "Stag C", "SR A", "Steph's B4", "Steph's C6", 
                                                       "SI A", "SI C", "SI D", "SI F", "Thicket-3", "TJ 2")]

########
#ANCOM-BC2 with genotypes as the group and time_point
########

# Subset to genotypes
rica_ancom_g = rica_ancom[, rica_ancom$a_Genotype %in% c("Acerv-1", "Acerv-2", "Acerv-3", "Cheetos A", "Cheetos D", "Cooper's",
                                                         "Elkhorn", "Kaufman-1", "Kelsey's 1", "KJS-74", "Marker 9",
                                                         "MB B", "Sandy", "Stag C", "SR A", "Steph's B4", "Steph's C6", 
                                                         "SI A", "SI C", "SI D", "SI F", "Thicket-3", "TJ 2")]


set.seed(123)

output = ancombc2(data = rica_ancom_g,
  assay_name = "counts",
  tax_level = "Family",
  fix_formula = "Time_Point + a_Genotype",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "a_Genotype",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = TRUE,
  iter_control = list(
    tol = 1e-2,
    max_iter = 20,
    verbose = TRUE ),
  em_control = list(
    tol = 1e-5,
    max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(
    fwer_ctrl_method = "holm",
    B = 100),
  trend_control = list(
    contrast = list(
      diag(22),  # Identity matrix as an example contrast
      matrix(c(rep(1, 22), rep(-1, 22)), nrow = 22, ncol = 22, byrow = TRUE),
      matrix(c(rep(-1, 22), rep(1, 22)), nrow = 22, ncol = 22, byrow = TRUE)),
    node = list(22, 22, 22),
    solver = "ECOS",
    B = 100))

#Structural zeros (Taxon presence/absence)
tab_zero = output$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")


########
#ANCOM-BC2 with genotypes as the group and spirocha y/n as a group
########

# Calculate variance for each taxon, assuming counts is a matrix or data frame
if (!is.null(rica_ancom$counts) && ncol(rica_ancom$counts) > 0) {
  taxa_variance <- apply(rica_ancom$counts, 1, var)
  
  # Filter out taxa with zero variance
  non_zero_variance_taxa <- names(taxa_variance[taxa_variance > 0])
  rica_ancom_filtered <- rica_ancom[non_zero_variance_taxa, ]
  
  # Print the filtered data for inspection
  print(rica_ancom_filtered)
} else {
  print("The counts matrix is empty or not correctly formatted.")
}


# Calculate variance for each taxon
taxa_variance <- apply(rica_ancom$counts, 1, var)

# Filter out taxa with zero variance
non_zero_variance_taxa <- names(taxa_variance[taxa_variance > 0])
rica_ancom_filtered <- rica_ancom[non_zero_variance_taxa, ]

# Print the filtered data for inspection
print(rica_ancom_filtered)

# Check the sample sizes for each group
table(rica_ancom_filtered$Spirochaetaceae_ANCOMpres)

# Remove groups with insufficient sample sizes (e.g., less than 5 samples)
sufficient_samples_groups <- names(table(rica_ancom_filtered$Spirochaetaceae_ANCOMpres)[table(rica_ancom_filtered$Spirochaetaceae_ANCOMpres) >= 5])
rica_ancom_filtered <- rica_ancom_filtered[, rica_ancom_filtered$Spirochaetaceae_ANCOMpres %in% sufficient_samples_groups]

# Print the filtered data for inspection
print(rica_ancom_filtered)

output2 <- ancombc2(
  data = rica_ancom_filtered,
  assay_name = "counts",
  tax_level = "ASVs",
  fix_formula = "Time_Point + Spirochaetaceae_ANCOMpres",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Spirochaetaceae_ANCOMpres",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = TRUE,
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
    B = 100
  ),
  trend_control = list(
    contrast = list(
      matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
      matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
      matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)
    ),
    node = list(2, 2, 1),
    solver = "ECOS",
    B = 10
  )
)

# Check the structure and summary of the output
str(output2)
summary(output2)

# Extracting LFCs for each taxon and factor combination
lfc_table <- output2$res

# Display the LFC table for inspection
print(lfc_table)

# Specific extraction of LFCs for 'Spirochaetaceae_ANCOMpres' categories
spiro_lfc <- lfc_table %>%
  select(taxon, contains("lfc_Spirochaetaceae_ANCOMpres"))

# Display the Spirochaetaceae specific LFCs
print(spiro_lfc)




set.seed(123)

# Subset to lean, overweight, and obese subjects
rica_ancom_s = rica_ancom[, rica_ancom$Spirochaetaceae_ANCOMpres %in% c("YES", "NO", "OTHER")]

# Note that by default, levels of a categorical variable in R are sorted 
# alphabetically. In this case, the reference level for `bmi` will be 
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
rica_ancom_s$Spirochaetaceae_ANCOMpres = factor(rica_ancom_s$Spirochaetaceae_ANCOMpres, levels = c("NO", "YES", "OTHER"))
print(rica_ancom_s)



output2 = ancombc2(data = rica_ancom,
                  assay_name = "counts",
                  tax_level = "ASVs",
                  fix_formula = "Time_Point + Spirochaetaceae_ANCOMpres",
                  rand_formula = NULL,
                  p_adj_method = "holm",
                  pseudo_sens = TRUE,
                  prv_cut = 0.10,
                  lib_cut = 1000,
                  s0_perc = 0.05,
                  group = "Spirochaetaceae_ANCOMpres",
                  struc_zero = TRUE,
                  neg_lb = TRUE,
                  alpha = 0.05,
                  n_cl = 2,
                  verbose = TRUE,
                  global = TRUE,  # Deactivated due to <3 categories
                  pairwise = TRUE,  # Deactivated due to <3 categories
                  dunnet = TRUE,  # Deactivated due to <3 categories
                  trend = TRUE,  # Deactivated due to <3 categories
                  iter_control = list(
                    tol = 1e-2,
                    max_iter = 20,
                    verbose = TRUE ),
                  em_control = list(
                    tol = 1e-5,
                    max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(
                    fwer_ctrl_method = "holm",
                    B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2, 1),
                                       solver = "ECOS",
                                       B = 10))

res_prim2 = output2$res

df_fig_TP2 = res_prim2 %>%
  dplyr::filter(diff_Spirochaetaceae_ANCOMpresYES == "TRUE") %>% 
  dplyr::arrange(desc(lfc_Spirochaetaceae_ANCOMpresYES)) %>%
  dplyr::mutate(direct = ifelse(lfc_Spirochaetaceae_ANCOMpresYES > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_Spirochaetaceae_ANCOMpresYES == "TRUE", "aquamarine3", "black"))


colnames(taxa)[1]<- "taxon"
df_fig_TP2 <- left_join(df_fig_TP2, taxa, by="taxon")
df_fig_TP2$taxon = factor(df_fig_TP2$taxon, levels = df_fig_TP2$taxon)

df_fig_TP2 %>% ggplot(aes(x = taxon, y = lfc_Spirochaetaceae_ANCOMpresYES, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin = lfc_Spirochaetaceae_ANCOMpresYES - se_Spirochaetaceae_ANCOMpresYES, 
                    ymax = lfc_Spirochaetaceae_ANCOMpresYES + se_Spirochaetaceae_ANCOMpresYES), 
                width = 0.2, position = position_dodge(0.05), color = "black")+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        legend.text = element_text(size = 15),
        title = element_text(size = 15))+
  labs(y="Log Fold Change", x="Taxa", title = "Log Fold Change from Spirochaetaceae (Presence/Absensense)")+
  coord_flip()+
  scale_x_discrete(labels= df_fig_TP2$ASVs)


##
#New Ancom Data
##
# Ensure Spirochaetaceae_ANCOMpres is a factor and has no missing values
rica_ancom$Spirochaetaceae_ANCOMpres <- factor(rica_ancom$Spirochaetaceae_ANCOMpres)
sum(is.na(rica_ancom$Spirochaetaceae_ANCOMpres))

# Run the ancombc2 function
output2 <- ancombc2(
  data = rica_ancom,
  assay_name = "counts",
  tax_level = "ASVs",
  fix_formula = "Time_Point + Spirochaetaceae_ANCOMpres",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Spirochaetaceae_ANCOMpres",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = FALSE,  # Set to FALSE due to <3 categories
  pairwise = FALSE,  # Set to FALSE due to <3 categories
  dunnet = FALSE,  # Set to FALSE due to <3 categories
  trend = FALSE,  # Set to FALSE due to <3 categories
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
    B = 100
  ),
  trend_control = list(
    contrast = list(
      matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
      matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE)
    ),
    node = list(2, 2, 1),
    solver = "ECOS",
    B = 10
  )
)

# Check the output
print(output2)

# Ensure Spirochaetaceae_ANCOMpres is a factor and has no missing values
rica_ancom$Spirochaetaceae_ANCOMpres <- factor(rica_ancom$Spirochaetaceae_ANCOMpres)
if (sum(is.na(rica_ancom$Spirochaetaceae_ANCOMpres)) > 0) {
  stop("There are missing values in Spirochaetaceae_ANCOMpres")
}

# Run the ancombc2 function
output2 <- ancombc2(
  data = rica_ancom,
  assay_name = "counts",
  tax_level = "ASVs",
  fix_formula = "Time_Point + Spirochaetaceae_ANCOMpres",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Spirochaetaceae_ANCOMpres",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = FALSE,  # Set to FALSE due to <3 categories
  pairwise = FALSE,  # Set to FALSE due to <3 categories
  dunnet = FALSE,  # Set to FALSE due to <3 categories
  trend = FALSE,  # Set to FALSE due to <3 categories
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
    B = 100
  ),
  trend_control = list(
    contrast = list(
      matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
      matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE)
    ),
    node = list(2, 2, 1),
    solver = "ECOS",
    B = 10
  )
)

# Check the structure and summary of the output
str(output2)
summary(output2)

# Extracting LFCs for each taxon and factor combination
lfc_table <- output2$res

# Display the LFC table for inspection
print(lfc_table)

# Specific extraction of LFCs for 'Spirochaetaceae_ANCOMpresYES'
spiro_lfc <- lfc_table[, c("taxon", "lfc_Spirochaetaceae_ANCOMpresYES")]

# Display the Spirochaetaceae specific LFCs
print(spiro_lfc)




########
#ANCOM-BC2 with just timepoint no grouping and using ASVs not by Family 
########

set.seed(123)

output3 = ancombc2(
  data = rica_ancom,
  assay_name = "counts",
  tax_level = "ASVs",
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


colnames(taxa)[1]<- "taxon"
df_fig_TP3 <- left_join(df_fig_TP3, taxa, by="taxon")
df_fig_TP3$taxon = factor(df_fig_TP3$taxon, levels = df_fig_TP3$taxon)

lfc_16 <- df_fig_TP3 %>% ggplot(aes(x = taxon, y = lfc_Time_PointT2, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin = lfc_Time_PointT2 - se_Time_PointT2, ymax = lfc_Time_PointT2 + se_Time_PointT2), 
                width = 0.2, position = position_dodge(0.05), color = "black")+
  scale_x_discrete(labels=c("Campylobacterales (751)", "Unknown Cyanobacteria (1235)", "SAR11_clade (1095)", "Gammaproteobacteria (980)",
                            "Pseudomonadaceae (1005)", "Salinimicrobium catena (274)", "Tenacibaculum (769)", "Bacteroidetes (1225)", "Endozoicomonas atrinae (740)",
                           "SAR86_clade (834)", "Unknown Cyanobacteria (1643)", "Paraburkholderia (683)", "Roseobacteraceae (531)", 
                           "Unknown Cyanobacteria (860)", "Unknown Cyanobacteria (546)", "Vibrio (341)", "Unknown Cyanobacteria (1024)", "Unknown Cyanobacteria (1010)",
                            "Unknown Cyanobacteria (1175)", "Unknown Cyanobacteria (1498)", "Unknown Cyanobacteria (729)", "Bacteroidetes (943)"))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        legend.text = element_text(size = 8),
        title = element_text(size = 8),
        axis.text.y = element_text( size=5))+
  labs(y="Log Fold Change", color="Log Fold Change", x="Taxa", title = "Log Fold Change from T1 to T2")+
  coord_flip()


##
#volcano plot for ancom 
##
res_prim3_lab <- res_prim3 %>%
  mutate(reg_status = case_when(lfc_Time_PointT2 < (-0.5)  ~ "Down Regulated",
                                lfc_Time_PointT2 > 0.5  ~ "Up Regulated",))



volc_16<- ggplot(res_prim3_lab, aes(x=lfc_Time_PointT2, y= -log10(p_Time_PointT2), col = reg_status))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed')+ 
  scale_color_manual(values = c( "#bb0c00","#00AFBB", "grey"), # to set the colours of our variable  
                     labels = c("Decrease", "Increase" , "Not significant"))+
  labs(title = "Volcano Plot for 16S data", 
       x = "Log Fold Change from T1 to T2", y = "-log10(p-value)",
       color="Abundance Status")+
  theme_classic2()+
  theme(axis.text.x = element_text( size=8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        legend.text = element_text(size = 8),
        title = element_text(size = 8),
        axis.text.y = element_text( size=8))


#combine 

Ancomn_Plots<-ggarrange(lfc_16, volc_16, font.label = list(size = 10, color = "black"), labels=c("A", "B"), common.legend = FALSE)

ggsave("/Volumes/RSMASFILES2/ThesisPlots/fig_6.png", width = 8, height = 4.5, units = "in",
       dpi=300)


