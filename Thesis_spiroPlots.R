library(tidyverse)

###
##Sterling Butler - NOAA - UM 
## 07/26/2024
## Spirochateceae enrichment with ACDC data
###



#import the meta data
meta_df<-read.csv("qiime2_05202024_16S_RICA_META.csv")

### filter out outlier
meta_df<- meta_df %>% filter(ratio_mean != "NA")

#get only one value for genotype
reduce_df <- meta_df %>% distinct(Genotype, .keep_all = TRUE)



#if you want to show what is significant and then what you want to compare
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
comparisons_material = list( c("YES", "NO"))


a <- ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=Growth_index))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15),
        legend.text = element_text(size = 5))+
  labs( 
       x = "Spirochaetaceae enriched", y = "Growth index ",
       color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)

a1 <- ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=Growth_6mon))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15),
        legend.text = element_text(size = 5))+
  labs( 
    x = "Spirochaetaceae enriched", y = "6 Month Growth Rate",
    color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)

a2 <- ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=Growth_Specific))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15),
        legend.text = element_text(size = 5))+
  labs( 
    x = "Spirochaetaceae enriched", y = "Specific Growth Rate",
    color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)



b <- ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=Heat_tolerance_CBASS))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15),
        legend.text = element_text(size = 5))+
  labs( 
    x = "Spirochaetaceae enriched", y = "Heat tolerance (CBASS)",
    color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)


c <-ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=Bleach_index))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15),
        legend.text = element_text(size = 5))+
  labs( 
    x = "Spirochaetaceae enriched", y = "Bleach index",
    color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)


d <-ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=Hcrit))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.title = element_text(size = 15),
        legend.text = element_text(size = 5))+
  labs( 
    x = "Spirochaetaceae enriched", y = "Hcrit ",
    color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)


e <-ggplot(reduce_df, aes(x=Spirochaetaceae_ANCOMpres, y=log(symbiont_density_hemocytometer)))+
  geom_boxplot(alpha = 2, aes(colour=Spirochaetaceae_ANCOMpres))+
  theme_classic2()+
  geom_jitter(size = .5)+
  scale_color_manual(values = c("#000000", "#7de661"))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.title = element_text(size = 12),
        legend.text = element_text(size = 5))+
  labs( 
    x = "Spirochaetaceae enriched", y = "symbiont density (hemocytometer)",
    color= "Spirochaetaceae enriched")  +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)



ggarrange(a, a1, a2, b, c, d, e, font.label = list(size = 10, color = "black"), labels=c("A", "B", "C", "D", "E", "F", "G"), 
          common.legend = TRUE)


ggsave("/Volumes/RSMASFILES2/ThesisPlots/spiro_acdc.png", width = 11, height = 8, units = "in",
       dpi=300)


#just growth variables boxplots

ggarrange(a, a1, a2, font.label = list(size = 10, color = "black"), labels=c("A", "B", "C", "D", "E", "F", "G"), 
          common.legend = TRUE)

ggsave("/Volumes/RSMASFILES2/ThesisPlots/spiro_acdc_growth.png", width = 11, height = 8, units = "in",
       dpi=300)



### wilcox test

wilcox.test(Growth_6mon~Spirochaetaceae_ANCOMpres , data = reduce_df,
                   exact = FALSE)

wilcox.test(symbiont_density_hemocytometer~Spirochaetaceae_ANCOMpres , data = reduce_df,
            exact = FALSE)

wilcox.test(Heat_tolerance_CBASS~Spirochaetaceae_ANCOMpres , data = reduce_df,
            exact = FALSE)

wilcox.test(Bleach_index~Spirochaetaceae_ANCOMpres , data = reduce_df,
            exact = FALSE)

wilcox.test(Hcrit~Spirochaetaceae_ANCOMpres , data = reduce_df,
            exact = FALSE)

#test for normalilty 

shapiro.test(meta_df$Growth_index)
shapiro.test(meta_df$Growth_6mon)
shapiro.test(meta_df$Growth_Specific)


### t-test

t.test(Growth_6mon~Spirochaetaceae_ANCOMpres , data = reduce_df)

t.test(Growth_Specific~Spirochaetaceae_ANCOMpres , data = reduce_df)

## anova models

anova_model <- aov(Growth_Specific ~ Spirochaetaceae_ANCOMpres, data = reduce_df)
summary(anova_model )

anova_model <- aov(Growth_index ~ Spirochaetaceae_ANCOMpres * Genotype, data = meta_df)
summary(anova_model )

anova_model <- aov(Growth_6mon~ Spirochaetaceae_ANCOMpres * Genotype, data = meta_df)
summary(anova_model )

library(car)
plot(anova_model, which = 3)


reduce_df <-reduce_df %>%
  mutate(Spirochaetaceae_ANCOMpres = if_else(Spirochaetaceae_ANCOMpres=="YES", 1, 0))

model <- glm(Spirochaetaceae_ANCOMpres~Growth_6mon, data= reduce_df, family = "binomial")
summary(model)

ggplot(reduce_df, aes(Growth_6mon, Spirochaetaceae_ANCOMpres)) +
  stat_smooth(method="glm", method.args = list(family = "binomial"), formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0))





