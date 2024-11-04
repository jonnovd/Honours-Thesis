#load packages
library(labdsv)
library(vegan)
library(ape)
library(phyloseq)
library(pgirmess)
library(dplyr)
library(gplots)
library(picante)
library(car)
library(ggplot2)
library(Hmisc)
library(RColorBrewer)
library(reshape2)
library(raincloudplots)
library(RAM)

setwd("")

#loading tables
taxa<-read.csv("Taxa_table.csv",header = T,row.names = 1, sep = ";")
meta<-read.csv("Meta_data.csv", header=T,row.names=1, sep = ";")
meta
library(ggplot2)

ncol(taxa) # number of ASVs
nrow(taxa) # number of samples
dim(taxa) # dimensions of data frame
row.names(taxa) # sample names
taxa[1:5,1:5] # show first 5 columns and first 5 rows
colnames(taxa) #colum/ASV names

## Filtering out underrepresented otus seen fewer than 5 times in the data
### These are uninformative and may be sequencing artifacts

taxa.sums<-colSums(taxa) # calculate total sums of every OTU
taxa.filtered<-taxa[, which(taxa.sums >5)] # get rid of OTUs seen less than 4 times

ncol(taxa.filtered) # how many OTUs do we have now?
ncol(taxa) # how many did we have before?

library(labdsv)
library(mgcv)

# Richness
n.otus<-rowSums(taxa.filtered >1) # number of ASVs in each sample
n.otus

###number of reads obtained per sample
rowSums(taxa.filtered)
rowSums(taxa)

summary(rowSums(taxa.filtered))
n.otus.rar<- rarefy(taxa.filtered, min(rowSums(taxa.filtered))) 
n.otus.rar

## ALPHA DIVERSITY
shann<-diversity(taxa.filtered) # calculate Shannon index of diversity
simp<-diversity(taxa.filtered, "invsimpson") # calculate Simpson index of diversity
div.taxa<-data.frame(n.otus, n.otus.rar, shann, simp)
div.taxa

# Perform rarefaction
# This ensures that all samples have the same sequencing depth (equal number of counts).
n.otus.rar <- rarefy(taxa.filtered, min(rowSums(taxa.filtered))) 

# Rarefy the counts to the same depth as determined by rarefaction
taxa.filtered.rar <- taxa.filtered
for (i in 1:nrow(taxa.filtered)) {
  taxa.filtered.rar[i, ] <- taxa.filtered[i, ] * (n.otus.rar[i] / rowSums(taxa.filtered)[i])
}

# Calculate alpha diversity using the rarefied dataset
shann.rar <- diversity(taxa.filtered.rar)  # Shannon index
simp.rar <- diversity(taxa.filtered.rar, "invsimpson")  # Simpson index

div.taxa <- data.frame(n.otus.rar, shann.rar, simp.rar)
div.taxa

#saving alpha diversity table
write.csv(div.taxa, "diversity.report.csv")

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(dplyr)

meta_tp <- meta[meta$Timepoint == "T6",]
meta_tp
simp_tp <- simp[names(simp) %in% rownames(meta_tp)]
simp_tp

meta_ab <- meta[meta$Treatment == "AB",]
meta_ab
#meta_ab <- meta_ab[meta_ab$Timepoint %in% c("T1", "T6"), ]
#meta_ab
simp 
simp.rar 

simp_ab <- simp[names(simp) %in% meta_ab$sample]
simp_ab

Simp.plot <- cbind.data.frame(simp_ab,meta_ab)
colnames(Simp.plot)
#Simp.plot <- cbind.data.frame(simp_tp,meta_tp)
p_value <- kruskal.test(simp_ab~meta_ab$Timepoint)$p.value
#p_value <- wilcox.test(simp_ab~meta_ab$Treatment)$p.value

#FRIEDMAN TEST
library(tidyverse)
library(ggpubr)
library(rstatix)
df_meta <- meta[meta$Treatment == "AB",]
df_meta
simpson <- simp[names(simp) %in% df_meta$sample]
simpson
Simp.plot <- cbind.data.frame(simpson,df_meta)

df_ftest <- cbind.data.frame(simpson,df_meta)
df_ftest$patient <- sapply(strsplit(df_meta$sample, "-"), `[`, 1)
df_ftest

res.fried <- df_ftest %>% friedman_test(simpson ~ Timepoint |patient)
res.fried

df_ftest %>% friedman_effsize(simpson ~ Timepoint |patient)

# Wilcoxon pairwise signed-rank tests for determining which comparisons are significanlty different:
pwc <- df_ftest %>%
  wilcox_test(simpson ~ Timepoint, paired = TRUE, p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% add_xy_position(x = "Timepoint")
pwc$Timepoint <- "T1"
pwc

ggplot(Simp.plot, aes(x=Timepoint, y=simpson, fill = Timepoint, shape = Timepoint)) + theme_bw() +
  geom_boxplot(outlier.color = "red", outlier.size = 3, outlier.shape = 8) +
  geom_point(position=position_jitterdodge()) + 
  scale_fill_manual(values=c("purple","turquoise", "lightgreen", "lightblue")) + 
  ylab("Inverse Simpson Diversity Index")+ 
  xlab(label = "Timepoint") +
  ggtitle(label = "Antibiotics") +
  scale_y_continuous(breaks = c(2.5, 5, 7.5, 10, 12.5)) +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.text.x=element_markdown(face="bold",size=10),
        axis.text.y=element_markdown(face="bold",size=14)) +
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none") 
  #stat_pvalue_manual(pwc, hide.ns = FALSE) +
  #labs(
   # subtitle = get_test_label(res.fried,  detailed = TRUE),
   # caption = get_pwc_label(pwc)
  #)

# CHANGE IN DIVERSITY STATISTICAL TEST

df_ab <- meta[meta$Treatment == "AB",]
df_ab
simpson_ab <- simp[names(simp) %in% df_ab$sample]
simpson_ab
df_simp_ab <- cbind.data.frame(simpson_ab,df_ab)
df_simp_ab <- rownames_to_column(df_simp_ab, var = "Sample")
df_simp_ab
df_simp_ab$Patient <- sapply(strsplit(df_simp_ab$Sample, "-"), `[`, 1)

df_simp_ab_t1 <- df_simp_ab[df_simp_ab$Timepoint == "T1", ]
df_simp_ab_t2 <- df_simp_ab[df_simp_ab$Timepoint == "T2", ]

df_simp_ab_t1 <- df_simp_ab_t1[order(df_simp_ab_t1$Patient), ]
df_simp_ab_t2 <- df_simp_ab_t2[order(df_simp_ab_t2$Patient), ]

# Calculate differences between T2 and T1 for each patient
df_simp_diff_ab <- data.frame(
  Patient = df_simp_ab_t1$Patient,
  T2minusT1 = df_simp_ab_t2$simpson_ab - df_simp_ab_t1$simpson_ab
)

df_simp_diff_ab

df_control <- meta[meta$Treatment == "CONTROL",]
df_control
simpson_control <- simp[names(simp) %in% df_control$sample]
simpson_control
df_simp_control <- cbind.data.frame(simpson_control,df_control)
df_simp_control <- rownames_to_column(df_simp_control, var = "Sample")
df_simp_control
df_simp_control$Patient <- sapply(strsplit(df_simp_control$Sample, "-"), `[`, 1)

df_simp_control_t1 <- df_simp_control[df_simp_control$Timepoint == "T1", ]
df_simp_control_t2 <- df_simp_control[df_simp_control$Timepoint == "T2", ]

df_simp_control_t1 <- df_simp_control_t1[order(df_simp_control_t1$Patient), ]
df_simp_control_t2 <- df_simp_control_t2[order(df_simp_control_t2$Patient), ]

# Calculate differences between T2 and T1 for each patient
df_simp_diff_control <- data.frame(
  Patient = df_simp_control_t1$Patient,
  T2minusT1 = df_simp_control_t2$simpson_control - df_simp_control_t1$simpson_control
)

df_simp_diff_control

wilcox_result <- wilcox.test(
  df_simp_diff_ab$T2minusT1, 
  df_simp_diff_control$T2minusT1,
  alternative = "less"  # Can be "two.sided", "greater", or "less"
)

# View the results
wilcox_result

# Load the ggplot2 library if not already loaded
library(ggplot2)

# Combine the difference data into one dataframe
df_simp_diff_ab$Treatment <- "Antibiotic"
df_simp_diff_control$Treatment <- "Control"
df_simp_diff_combined <- rbind(df_simp_diff_ab, df_simp_diff_control)

# Boxplot with ggplot2
p <- ggplot(df_simp_diff_combined, aes(x = Treatment, y = T2minusT1, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Change in Diversity from Baseline to Treatment Cessation",
       y = "Change in Inverse Simpson Diversity Index (T2 - T1)") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_markdown(face = "bold", size = 10),
    axis.text.y = element_markdown(face = "bold", size = 14),
    axis.title.x = element_markdown(face = "bold", size = 14),
    axis.title.y = element_markdown(face = "bold", size = 14)
  )

# Add significance annotation based on Wilcoxon test result
p + annotate("text", x = 1.5, y = max(df_simp_diff_combined$T2minusT1) * 1.1,
             label = paste("Wilcoxon p-value:", round(wilcox_result$p.value, 3)),
             size = 4)


# BETA DIVERSITY ANALYSIS
#PCOA USING GGPLOT

ggplot(data = as.data.frame(pcoa.taxa$vectors[,1:2]), aes(y=Axis.2,x=Axis.1,color = meta_ab$Timepoint),cex=3,pch=21, add=TRUE) +
  #geom_point(size = 5)+ 
  geom_point(aes(y = Axis.2, x = Axis.1, fill = meta_tp$Treatment), 
             size = 5, shape = 21, color = "black", stroke = 0.3, 
             position = position_jitter(width = 0.02, height = 0.02)) +
  #scale_fill_manual(values = c("AB" = "lightcoral", "CONTROL" = "lightblue")) +
  ggtitle(label = "Bray-Curtis B Diversity at T1 (baseline)")+
  xlab(label = "Percent Variation at 64%")+
  ylab(label = "Percent Variation at 15%")+
  xlim(-0.4, 0.7) +
  ylim(-0.5, 0.5) +
  labs(fill = "Treatment")+
  #stat_ellipse()+
  theme_bw()+
  #scale_color_manual(values=c("orange","green","red","blue","lavender","purple"))+
  #scale_fill_manual(values=c("orange","green","red","blue","lavender","purple")) + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty =2)+
  geom_hline(yintercept = 0, lty =2) + 
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.x=element_markdown(face="bold",size=8),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(face="bold",size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14))

# SETTINGS FOR T2
meta_tp <- meta[meta$Timepoint == "T2",]
meta_tp
rownames(meta_tp)
rownames(taxa.filtered)
taxa.filtered.tp <- taxa.filtered[rownames(taxa.filtered) %in% meta_tp$sample, ]
dim(taxa.filtered.tp)
taxa.filtered.tp.relab<-decostand(taxa.filtered.tp, "total")*100
dim(taxa.filtered.tp.relab)

library(ape)
taxa.bray.relab<-vegdist(taxa.filtered.tp.relab, method = "bray")
dim(taxa.bray.relab)
pcoa.taxa<-pcoa(taxa.bray.relab)
pcoa.taxa$values
biplot(pcoa.taxa)

#PCOA USING GGPLOT
ggplot(data = as.data.frame(pcoa.taxa$vectors[,1:2]), aes(y=Axis.2,x=Axis.1,color = meta_tp$Treatment),cex=3,pch=21, add=TRUE) +
  #geom_point(size = 5)+ 
  geom_point(aes(y = Axis.2, x = Axis.1, fill = meta_tp$Treatment), 
             size = 5, shape = 21, color = "black", stroke = 0.3, 
             position = position_jitter(width = 0.02, height = 0.02)) +
  #scale_fill_manual(values = c("AB" = "lightcoral", "CONTROL" = "lightblue")) +
  ggtitle(label = "Bray-Curtis B Diversity at T2 (treatment cessation)")+
  xlab(label = "Percent Variation at 40%")+
  ylab(label = "Percent Variation at 21%")+
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.7) +
  labs(fill = "Treatment")+
  #stat_ellipse()+
  theme_bw()+
  #scale_color_manual(values=c("orange","green","red","blue","lavender","purple"))+
  #scale_fill_manual(values=c("orange","green","red","blue","lavender","purple")) + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty =2)+
  geom_hline(yintercept = 0, lty =2) + 
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.x=element_markdown(face="bold",size=8),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(face="bold",size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14))

# SETTINGS FOR T6
meta_tp <- meta[meta$Timepoint == "T6",]
meta_tp
rownames(meta_tp)
rownames(taxa.filtered)
taxa.filtered.tp <- taxa.filtered[rownames(taxa.filtered) %in% meta_tp$sample, ]
dim(taxa.filtered.tp)
taxa.filtered.tp.relab<-decostand(taxa.filtered.tp, "total")*100
dim(taxa.filtered.tp.relab)

library(ape)
taxa.bray.relab<-vegdist(taxa.filtered.tp.relab, method = "bray")
dim(taxa.bray.relab)
pcoa.taxa<-pcoa(taxa.bray.relab)
pcoa.taxa$values
biplot(pcoa.taxa)

#PCOA USING GGPLOT
ggplot(data = as.data.frame(pcoa.taxa$vectors[,1:2]), aes(y=Axis.2,x=Axis.1,color = meta_tp$Treatment),cex=3,pch=21, add=TRUE) +
  #geom_point(size = 5)+ 
  geom_point(aes(y = Axis.2, x = Axis.1, fill = meta_tp$Treatment), 
             size = 5, shape = 21, color = "black", stroke = 0.3, 
             position = position_jitter(width = 0.02, height = 0.02)) +
  #scale_fill_manual(values = c("AB" = "lightcoral", "CONTROL" = "lightblue")) +
  ggtitle(label = "Bray-Curtis B Diversity at T6 (6 months corrected age)")+
  xlab(label = "Percent Variation at 40%")+
  ylab(label = "Percent Variation at 23%")+
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  labs(fill = "Treatment")+
  #stat_ellipse()+
  theme_bw()+
  #scale_color_manual(values=c("orange","green","red","blue","lavender","purple"))+
  #scale_fill_manual(values=c("orange","green","red","blue","lavender","purple")) + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty =2)+
  geom_hline(yintercept = 0, lty =2) + 
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.x=element_markdown(face="bold",size=8),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(face="bold",size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14))


# SETTINGS FOR Antibiotic group
meta_ab <- meta[meta$Treatment == "AB",]
meta_ab <- rownames_to_column(meta_ab, var = "sample")
meta_ab
rownames(meta_ab)
rownames(taxa.filtered)
taxa.filtered.ab <- taxa.filtered[rownames(taxa.filtered) %in% meta_ab$sample, ]
dim(taxa.filtered.ab)
taxa.filtered.ab.relab<-decostand(taxa.filtered.ab, "total")*100
dim(taxa.filtered.ab.relab)
taxa.bray.relab<-vegdist(taxa.filtered.ab.relab, method = "bray")
dim(taxa.bray.relab)
pcoa.taxa<-pcoa(taxa.bray.relab)
pcoa.taxa$values
biplot(pcoa.taxa)

#PCOA USING GGPLOT
ggplot(data = as.data.frame(pcoa.taxa$vectors[,1:2]), aes(y=Axis.2,x=Axis.1,color = meta_ab$Timepoint),cex=3,pch=21, add=TRUE) +
  #geom_point(size = 5)+ 
  geom_point(aes(y = Axis.2, x = Axis.1, fill = meta_ab$Timepoint), 
             size = 5, shape = 21, color = "black", stroke = 0.3, 
             position = position_jitter(width = 0.03, height = 0.02)) +
  #scale_fill_manual(values = c("AB" = "lightcoral", "CONTROL" = "lightblue")) +
  ggtitle(label = "Bray-Curtis B Diversity - Antibiotics")+
  xlab(label = "Percent Variation at 28%")+
  ylab(label = "Percent Variation at 13%")+
  xlim(-0.75, 1.0) +
  ylim(-0.75, 0.5) +
  labs(fill = "Timepoint")+
  #stat_ellipse()+
  theme_bw()+
  #scale_color_manual(values=c("orange","green","red","blue","lavender","purple"))+
  #scale_fill_manual(values=c("orange","green","red","blue","lavender","purple")) + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty =2)+
  geom_hline(yintercept = 0, lty =2) + 
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.x=element_markdown(face="bold",size=8),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(face="bold",size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14))

#Relative Abundance
dim(taxa.filtered)
colnames(taxa.filtered)

taxa.filtered.relab<-decostand(as.data.frame(taxa.filtered), "total")*100
taxa.filtered.relab.2 <- sweep(t(taxa.filtered.relab),2,colSums(t(taxa.filtered.relab)),'/')
taxa.filtered.relab.2 <- taxa.filtered.relab.2[rowMeans(taxa.filtered.relab.2) >= .01,]

taxa.filtered.relab.2 <- taxa.filtered.relab.2[order(rowMeans(taxa.filtered.relab.2),decreasing = F),]

plot <- as.data.frame(t(taxa.filtered.relab.2))

plot <- rownames_to_column(plot, var = "sample")
colnames(plot)

library(reshape2)
plot <- melt(plot, id = "sample", variable.name = "Species")
colnames(meta)
plot <- merge(plot, meta, by = "sample")

plot <- plot %>%
  group_by(Treatment, Timepoint, Species) %>%
  dplyr::summarise(mean_abundance = mean(value, na.rm = TRUE))

# Calculate the "Other" category for species with low abundance
other <- plot %>%
  group_by(Treatment, Timepoint) %>%
  dplyr::summarise(mean_abundance = 1 - sum(mean_abundance)) %>%
  mutate(Species = "Other")

# Combine the main species data with the "Other" category
plot <- as.data.frame(plot)
other <- as.data.frame(other)
plot <- rbind(plot, other)

# Reorder species
tax_ord_factor <- as.character(rownames(taxa.filtered.relab.2))
plot$Species <- factor(plot$Species, levels = tax_ord_factor)
plot$Species <- gsub('\\.', ' ', plot$Species)
plot$Species[is.na(plot$Species)] <- "Other"

plot$Species <- reorder(plot$Species, plot$mean_abundance)
plot$Species <- factor(plot$Species, levels=rev(levels(plot$Species)))

# Create the bar plot with ggplot
ggplot(plot, aes(fill = Species, y = mean_abundance, x = Timepoint)) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity", width = 0.5) +
  scale_fill_manual(values = c(
  "#1f78b4", # S marcescens
  "grey",
  "darkgreen",
  "lightgreen",
  "blue",
  "#ff7f00",
  "red", 
  "pink",
  "#ff33a3",
  "purple",
  "yellow2",
  "turquoise",
  "#66a61e",
  "coral3",
  "lightblue",
  "#ffcc00",
  "#d95f02",
  "#7570b3",
  "#e7298a",
  "#ffff99",
  "lightgreen",
  "lightblue",
  "blue",
  "purple",
  "pink",
  "red",
  "orange")) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 1, byrow = FALSE, face = "italic")) +
  xlab("Timepoint") +
  ylab("Mean Relative Abundance (%)") +
  #labs(title = "Mean Microbial Relative Abundance", fill = "Species") +
  theme(
    axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_markdown(face = "bold", size = 13),
    legend.text = element_text(size = 14, face = "italic"),
    legend.title = element_text(face = "bold", size = 15),
    axis.title.x = element_markdown(face = "bold", size = 16, margin = margin(t = 10)),
    axis.title.y = element_markdown(face = "bold", size = 16, margin = margin(r = 10)),
    plot.title = element_text(face = "bold", size = 20),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  facet_grid(~Treatment, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 10))

# ATTEMPT AT LINEAR MIXED EFFECTS MODEL
# install.packages("lme4")
# install.packages("lmerTest")  # For p-values on fixed effects

# library(lme4)
# library(lmerTest)
# install.packages("tibble")
# library(tibble)

# div.taxa
# # Assuming your data frame is named 'div.taxa'
# div.taxa <- rownames_to_column(div.taxa, var = "sample")

# # View the updated data frame
# head(div.taxa)
# div.taxa
# meta
# meta_ab <- meta[meta$Treatment == "AB",]

# # Filter data for the antibiotic treatment group
# data_antibiotic <- subset(div.taxa, sample %in% meta_ab$sample)
# data_antibiotic
# data_antibiotic <- merge(data_antibiotic, meta_ab[, c("sample", "Timepoint")], by = "sample")
# data_antibiotic
# data_antibiotic$sample_id <- sub("-.*", "", data_antibiotic$sample)
# data_antibiotic

# # Linear mixed model for Simpson index with Timepoint as a fixed effect
# lmm_simp <- lme4::lmer(simp ~ Timepoint + (1|sample_id), data = data_antibiotic)
# summary(lmm_simp)
# library(emmeans)

# # Post-hoc test for pairwise comparisons between time points
# posthoc_simp <- emmeans(lmm_simp, pairwise ~ Timepoint)
# summary(posthoc_simp)

# library(ggplot2)

# emmeans_values <- as.data.frame(emmeans(lmm_simp, ~ Timepoint))
# # Boxplot of Simpson diversity across timepoints
# ggplot(data_antibiotic, aes(x = Timepoint, y = simp, fill = Timepoint)) +
#   geom_boxplot() +
#   geom_point(data = emmeans_values, aes(x = Timepoint, y = emmean), 
#              size = 3, shape = 23, fill = "red") +  # Adds emmeans
#   theme_bw() + 
#   labs(x = "Timepoint", y = "Simpson Diversity Index", 
#        title = "Simpson Diversity Across Timepoints") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(face = "bold", size = 12),
#         axis.text.y = element_text(face = "bold", size = 12),
#         axis.title = element_text(face = "bold", size = 14))

# HEATMAP
library(reshape2)
library(pheatmap)

# Function to generate a heatmap for a specific treatment group, sorted by timepoints
generate_heatmap_by_treatment <- function(data, treatment_group) {
  # Reshape the data: Convert from long to wide format
  heatmap_data_long <- melt(data, id.vars = c("sample", "Treatment", "Timepoint"),
                            variable.name = "Species", value.name = "abundance")
  
  # Cast the data: Combine Sample and Timepoint in rows, species in columns
  heatmap_matrix_filtered <- dcast(heatmap_data_long, paste(sample, Timepoint, sep = "_") ~ Species, value.var = "abundance")
  
  # Set rownames to a combination of Sample and Timepoint
  rownames(heatmap_matrix_filtered) <- heatmap_matrix_filtered[, 1]  # First column contains Sample_Timepoint
  heatmap_matrix_filtered <- heatmap_matrix_filtered[, -1]  # Remove the combined Sample_Timepoint column
  
  species_order <- c("Serratia marcescens", "Moraxella catarrhalis", "Moraxella nonliquefaciens", 
                     "Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus epidermidis", 
                     "Cutibacterium acnes", "Enterobacter hormaechei", "Streptococcus mitis",
                     "Staphylococcus aureus", "Streptococcus pneumoniae", "Streptococcus oralis",
                     "Dolosigranulum pigrum", "Gemella haemolysans", "Rothia mucilaginosa")
  
  heatmap_matrix_filtered <- heatmap_matrix_filtered[, species_order, drop = FALSE]
 
  # Extract sample names and timepoints
  samples <- sapply(strsplit(rownames(heatmap_matrix_filtered), "-"), `[`, 1)
  timepoints <- sapply(strsplit(rownames(heatmap_matrix_filtered), "_"), `[`, 2)
  
  # Create a dataframe for annotation that separates samples
  annotation_df <- data.frame(Timepoint = timepoints)
  rownames(annotation_df) <- rownames(heatmap_matrix_filtered)
  
  # Define colors for each timepoint (this will visually separate timepoints)
  timepoint_colors <- list(Timepoint = c("T1" = "white", "T2" = "grey80", "T4" = "grey60", "T6" = "grey40"))
  # Calculate relative abundance per sample (column-wise normalization)
  heatmap_matrix_relative <- sweep(heatmap_matrix_filtered, 1, rowSums(heatmap_matrix_filtered), "/")
  
  # Generate the heatmap without row scaling
  pheatmap(
    heatmap_matrix_relative,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = annotation_df,
    annotation_colors = timepoint_colors,
    show_rownames = FALSE,     # Hide sample names
    show_colnames = TRUE,     # Show species names
    fontsize_col = 12,
    color = colorRampPalette(c("grey97", "coral", "firebrick3"))(50),  # Define color scale
    #main = sprintf("Relative Abundance Heatmap of %s Group", treatment_group),
    gaps_row = cumsum(table(samples)),  # Create gaps between different samples
    border_color = "black",
  )
   
}
taxa.filtered
taxa.filtered.relab<-decostand(as.data.frame(taxa.filtered), "total")*100
taxa.filtered.relab.2 <- sweep(t(taxa.filtered.relab),2,colSums(t(taxa.filtered.relab)),'/')
taxa.filtered.relab.2 <- taxa.filtered.relab.2[rowMeans(taxa.filtered.relab.2) >= .002,]

taxa.filtered.relab.2 <- taxa.filtered.relab.2[order(rowMeans(taxa.filtered.relab.2),decreasing = F),]

num_species <- nrow(taxa.filtered.relab.2)
num_species

# Set the number of top species to plot, adjusting if there are fewer than 30
top_n <- min(15, num_species)

# Filter for top N species based on row means
top_species <- rownames(taxa.filtered.relab.2)[(num_species - top_n + 1):num_species]

# Filter the dataset to only include the top N species
taxa_top_n <- taxa.filtered.relab.2[rownames(taxa.filtered.relab.2) %in% top_species, ]

# Create a data frame for heatmap
heatmap_data <- as.data.frame(t(taxa_top_n))
heatmap_data <- rownames_to_column(heatmap_data, var = "sample")
colnames(heatmap_data) <- gsub("\\.", " ", colnames(heatmap_data))
colnames(heatmap_data)
heatmap_data <- merge(heatmap_data, meta[, c("sample", "Treatment", "Timepoint")], by = "sample")
colnames(heatmap_data)

dim(heatmap_data)
# Subset the data by Treatment group (Control and Antibiotics)
control_data <- subset(heatmap_data, Treatment == "CONTROL")
antibiotics_data <- subset(heatmap_data, Treatment == "AB")

# Generate heatmap for Control group, sorted by timepoints
generate_heatmap_by_treatment(control_data, "Control")

# Generate heatmap for Antibiotics group, sorted by timepoints
generate_heatmap_by_treatment(antibiotics_data, "Antibiotics")


# Antimicrobial Resistance Analysis (ARG) using data from ABRicate
amr_table <- read.csv("summary.tab",header = T,row.names = 1, sep = "\t", check.names=FALSE)
amr_table <- as.data.frame(amr_table)
amr_table$Patient <- sapply(strsplit(sapply(strsplit(rownames(amr_table), "-"), `[`, 1), "/"), `[`, 2)
amr_table
# Load the libraries
library(gt)
library(webshot2)

# Convert your data frame to a gt table
amr_table_gt <- gt(amr_table)
#selected_columns <- amr_table[, c("Patient", "aac(3)-XI", "aac(6')_Serra", "ant(9)-Ib", "aph(3')-Ia", "blaEC", "blaSST-1", "blaI_of_Z", "blaR1", "blaSST-1", "blaZ"), drop = FALSE]
#selected_columns <- amr_table[, c("Patient", "catA7", "dfrC", "erm(B)", "erm(C)", "erm(X)", "fosA5", "fusB", "mecR1", "mef(A)", "mph(C)", "msr(A)", "msr(D)"), drop = FALSE]
selected_columns <- amr_table[, c("Patient", "oqxB5", "qnrS1", "tet(41)", "tet(K)", "tet(M)"), drop = FALSE]
amr_table_gt <- gt(selected_columns)

amr_table_gt <- amr_table_gt %>%
  fmt_number(
    columns = everything(),
    decimals = 2
  ) %>%
  cols_move_to_start(
    columns = "Patient"
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "grey",
      weight = px(1)
    ),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(
      align = "center"
    ),
    locations = cells_body(columns = everything())
  ) %>%
  tab_options(
    table.width = px(800),
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  ) %>%
  cols_width(
    everything() ~ px(80)
  )

amr_table_gt
