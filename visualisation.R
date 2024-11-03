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

setwd("~/Documents/Honours/thesis/pipeline/public-data/preterm/diversity/")

#loading tables
taxa<-read.csv("Taxa_table.csv",header = T,row.names = 1, sep = ";")
meta<-read.csv("Meta_data.csv", header=T,row.names=1, sep = ";")
meta

#cleaning up

# grouped bar plot using ggplot

library(ggplot2)


ncol(taxa) # number of ASVs
nrow(taxa) # number of samples
dim(taxa) # dimensions of data frame
row.names(taxa) # sample names
taxa[1:5,1:5] # show first 5 columns and first 5 rows
colnames(taxa) #colum/ASV names

##Filtering - underrepresented otus, seen fewer than 4 times in the dataset 
###are uninformative and may be sequencing artifacts

taxa.sums<-colSums(taxa) # calculate total sums of every OTU
taxa.filtered<-taxa[, which(taxa.sums >5)] # get rid of OTUs seen less than 4 times

ncol(taxa.filtered) # how many OTUs do we have now?
ncol(taxa) # how many did we have before?


###Remove asvs present in only 1 samples or fewer
library(labdsv)
library(mgcv)
#taxa.filtered<-dropspc(taxa,2)
ncol(taxa.filtered)

# Let’s calculate richness in each group….
n.otus<-rowSums(taxa.filtered >1) # number of ASVs in each sample
n.otus

###number of reads obtained per sample
rowSums(taxa.filtered)
rowSums(taxa)


##calculates number of OTUs observed controlling for the minimum number of reads obtained in one given sample
summary(rowSums(taxa.filtered))


n.otus.rar<- rarefy(taxa.filtered, min(rowSums(taxa.filtered))) 
n.otus.rar


###Calculating Alpha Diversity
shann<-diversity(taxa.filtered) # calculate Shannon index of diversity
simp<-diversity(taxa.filtered, "invsimpson") # calculate Simpson index of diversity
div.taxa<-data.frame(n.otus, n.otus.rar, shann, simp)
div.taxa

# Step 1: Perform rarefaction
# This ensures that all samples have the same sequencing depth (equal number of counts).
n.otus.rar <- rarefy(taxa.filtered, min(rowSums(taxa.filtered))) 

# Step 2: Rarefy the counts to the same depth as determined by rarefaction
taxa.filtered.rar <- taxa.filtered
for (i in 1:nrow(taxa.filtered)) {
  taxa.filtered.rar[i, ] <- taxa.filtered[i, ] * (n.otus.rar[i] / rowSums(taxa.filtered)[i])
}

# Step 3: Calculate alpha diversity using the rarefied dataset
shann.rar <- diversity(taxa.filtered.rar)  # Shannon index
simp.rar <- diversity(taxa.filtered.rar, "invsimpson")  # Simpson index

# Step 4: Create a dataframe to store the results
div.taxa <- data.frame(n.otus.rar, shann.rar, simp.rar)
div.taxa


#saving alpha diversity table
write.csv(div.taxa, "diversity.report.csv")

meta

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(dplyr)
#meta_t1 <- meta %>% filter(Timepoint == "T1")
#samples_t1 <- rownames(meta_t1)
#samples_t1
#df_filtered <- div.taxa %>% filter(rownames(df) %in% samples_t1)
# Filter the shann dataframe to include only columns with "t1" in their names

meta_t1 <- meta[meta$Timepoint == "T1",]
meta_t1
shann_t1 <- shann[names(shann) %in% rownames(meta_t1)]
shann_t1

meta_ab <- meta[meta$Treatment == "AB",]
meta_ab
meta_ab <- rownames_to_column(meta_ab, var = "sample")
meta_ab
#meta_ab <- meta_ab[meta_ab$Timepoint %in% c("T1", "T2"), ]
#meta_ab
shann_ab <- shann[names(shann) %in% rownames(meta_ab)]
shann_ab

#Shann.plot <- cbind.data.frame(shann_t1,meta_t1)
Shann.plot <- cbind.data.frame(shann_ab,meta_ab)
p_value <- kruskal.test(shann_ab, meta_ab$Timepoint)$p.value

ggplot(Shann.plot, aes(x=Timepoint, y=shann_ab, fill = Timepoint, shape = Timepoint)) + theme_bw() +geom_boxplot() +
  geom_boxplot(outlier.color = "red", outlier.size = 3, outlier.shape = 8) +
  geom_point(position=position_jitterdodge()) + 
  #scale_color_manual(values=c("purple","turquoise", "lightgreen", "lightblue"))+
  #scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=c("purple","turquoise", "lightgreen", "lightblue")) + 
  ylab("Shannon Diversity")+ 
  xlab(label = "Timepoint") +
  ggtitle(label = "Control") +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.text.x=element_markdown(face="bold",size=10),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14)) +
  annotate("text", x = 1.5, y = max(Shann.plot$shann) + 0.3, 
           label = paste("p-value =", round(p_value, 4)), size = 4, hjust = 0, color="black")

#statistics
#kruskal.test(shann,meta$Timepoint)

#wilcox.test(shann~meta$Timepoint)
#wilcox.test(shann_t1~meta_t1$Treatment)



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




# Unweighted B diversity
taxa.bray.unw<-vegdist(decostand(taxa.filtered.ab.relab, "pa"), method="bray")
pcoa.taxa.unw<-pcoa(taxa.bray.unw)
pcoa.taxa.unw$values

#create a pcoa plot
biplot(pcoa.taxa.unw)
#unweighted
ggplot(data = as.data.frame(pcoa.taxa.unw$vectors[,1:2]), aes(y=Axis.2,x=Axis.1,color = meta_ab$Timepoint),cex=3,pch=21, add=TRUE) +
  geom_point(size = 4)+ 
  ggtitle(label = "UNW:Beta Diversity/Index:Bray-Curtis")+
  xlab(label = "Percent Variation at 15%")+
  ylab(label = "Percent Variation at 10%")+
  labs(color = "Treatment")+
  stat_ellipse()+
  theme_bw()+
  #scale_color_manual(values=c("orange","green","red","blue","lavender","purple"))+
  #scale_fill_manual(values=c("orange","green","red","blue","lavender","purple")) + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty =2)+
  geom_hline(yintercept = 0, lty =2) + 
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.text.x=element_markdown(face="bold",size=10),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(face="bold",size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14)) 


#Relative Abundance
dim(taxa.filtered)
colnames(taxa.filtered)
#taxa.filtered.no.Sm <- taxa.filtered %>% select(-Serratia.marcescens)

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







# stacked SPECIES Mycobacterium
##Phylum Stacked Bar####
#otu.phy<-t(read.table("phylum_taxa_table.tsv", header=T, row.names = 1,sep = "\t"))
#phy.filtered<-dropspc(otu.phy, 5)
taxa.filtered.relab<-decostand(as.data.frame(taxa.filtered), "total")*100
taxa.filtered.relab.2 <- sweep(t(taxa.filtered.relab),2,colSums(t(taxa.filtered.relab)),'/')
taxa.filtered.relab.2 <- taxa.filtered.relab.2[rowMeans(taxa.filtered.relab.2) >= .003,]

taxa.filtered.relab.2 <- taxa.filtered.relab.2[order(rowMeans(taxa.filtered.relab.2),decreasing = F),]

plot <- as.data.frame(t(taxa.filtered.relab.2))
plot <- rownames_to_column(plot, var = "Sample")

library(reshape2)
plot <- melt(plot, id = "Sample", variable.name = "Species")
plot <- plot %>% group_by(Sample, Species)
plot <- plot %>% dplyr::summarise(newvalue = sum(value))

other <- plot %>% group_by(Sample) %>% dplyr::summarise(newvalue = 1 - sum(newvalue))
other$Species <- "Other"
other <- other %>% dplyr::select(Sample, Species, newvalue)

plot <- as.data.frame(plot)
other <- as.data.frame(other)
plot <- as.data.frame(rbind(plot,other))

meta$Sample <- rownames(meta)
plot <- merge(plot, meta, by = "Sample")

tax_ord_factor<-as.character(rownames(taxa.filtered.relab.2))

plot$Species <- as.factor(plot$Species)
plot$Species <- factor(plot$Species, levels = tax_ord_factor)
plot <- plot[complete.cases(plot$Species),]

plot$Species <- gsub('\\.', ' ', plot$Species)


plot$Species <- reorder(plot$Species, plot$newvalue)
plot$Species <- factor(plot$Species, levels=rev(levels(plot$Species)))


#for alphabetically ordering the facetgrp
# Assuming `plot` is your data frame
#plot$Group1 <- factor(plot$Group1, levels=rev(levels(plot$Group1)))
#plot$Group1<-reorder(plot$Group1, plot$Group1)
#plot$Group1 <- factor(plot$Group1, levels = sort(unique(plot$Group1)))
#plot$Group2 <- factor(plot$Group2, levels = sort(unique(plot$Group2)))
#plot$Group1_Group2 <- interaction(plot$Group1, plot$Group2)





ggplot(plot, aes(fill=Species, y=newvalue, x=Sample)) +
  geom_bar(position= position_fill(reverse = TRUE), stat = "identity",  width = 0.5)+
  #scale_fill_manual(values = c(
    #"#1f78b4", # Proteobacteria
    #"#33a02c", # Bacteroidetes
    #"#e31a1c", # Actinobacteria
    #"#ff7f00", # Firmicutes
    #"#6a3d9a", # Cyanobacteria
    #"#b15928", # Verrucomicrobia
    #"#ffff99", # Chloroflexi
    #"#a6cee3", # Planctomycetes
    #"#fb9a99", # Nitrospirae
    #"#cab2d6", # Candidatus Gracilibacteria
    #"#fdbf6f", # Gemmatimonadetes
    #"#b2df8a", # Fusobacteria
    #"#ff33a3", # Euryarchaeota
    #"#c994c7", # Candidatus Saccharibacteria
    #"#bc80bd", # Acidobacteria
    #"#8dd3c7", # Mint Green
    #"#ffcc00", # Bright Yellow
    #"#d95f02", # Dark Orange
    #"#7570b3", # Purple
    #"#e7298a", # Hot Pink
    #"#66a61e" )) +
  theme_bw()+
  guides(fill=guide_legend(ncol=1, byrow=FALSE, face="italic"))+
  xlab("")+
  theme(legend.title = element_text("levels"),
        legend.key.height = unit(0.1, "mm"),
        legend.key.size = unit(10, "pt")) +
  theme_bw() +
  labs(title = "Bacterial pathogens", x= " group")+
  ylab(expression("Mean relative abundance(%)")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12), # Adjust angle and other properties
        #theme(axis.text.x=element_markdown(face="bold",size=12),
        axis.text.y=element_markdown(face="bold",size=20),
        legend.text=element_text(size=14, face="italic"),
        legend.title=element_text(face="bold",size=20)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=20),
        axis.title.y=element_markdown(face="bold",size=20)) +
  theme(plot.title = element_text(face="bold",size=20)) +
  facet_grid(~Treatment, scales = "free_x", space = "free_x")+
  theme(strip.text.x = element_text(size = 10))



# Trying Linear Mixed Effects model

install.packages("lme4")
install.packages("lmerTest")  # For p-values on fixed effects

library(lme4)
library(lmerTest)
install.packages("tibble")
library(tibble)

div.taxa
# Assuming your data frame is named 'div.taxa'
div.taxa <- rownames_to_column(div.taxa, var = "sample")

# View the updated data frame
head(div.taxa)
div.taxa
meta
meta_ab <- meta[meta$Treatment == "AB",]
meta_ab
#meta_ab <- rownames_to_column(meta_ab, var = "sample")
meta_ab

# Filter data for the antibiotic treatment group
data_antibiotic <- subset(div.taxa, sample %in% meta_ab$sample)
data_antibiotic
data_antibiotic <- merge(data_antibiotic, meta_ab[, c("sample", "Timepoint")], by = "sample")
data_antibiotic
data_antibiotic$sample_id <- sub("-.*", "", data_antibiotic$sample)
data_antibiotic

# Linear mixed model for Simpson index with Timepoint as a fixed effect
lmm_simp <- lme4::lmer(simp ~ Timepoint + (1|sample_id), data = data_antibiotic)

# Summary of the model
summary(lmm_simp)

#install.packages("emmeans")
library(emmeans)

# Post-hoc test for pairwise comparisons between time points
posthoc_simp <- emmeans(lmm_simp, pairwise ~ Timepoint)
summary(posthoc_simp)

library(ggplot2)

emmeans_values <- as.data.frame(emmeans(lmm_simp, ~ Timepoint))
# Boxplot of Simpson diversity across timepoints
ggplot(data_antibiotic, aes(x = Timepoint, y = simp, fill = Timepoint)) +
  geom_boxplot() +
  geom_point(data = emmeans_values, aes(x = Timepoint, y = emmean), 
             size = 3, shape = 23, fill = "red") +  # Adds emmeans
  theme_bw() + 
  labs(x = "Timepoint", y = "Simpson Diversity Index", 
       title = "Simpson Diversity Across Timepoints") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 14))



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
  
  #ordered_indices <- order(rownames(heatmap_matrix_filtered))
  #heatmap_matrix_filtered <- heatmap_matrix_filtered[ordered_indices, ]
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
  #annotation_df <- data.frame(Patient = samples, Timepoint = timepoints)
  #rownames(annotation_df) <- rownames(heatmap_matrix_filtered)
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
meta
#meta$Sample <- rownames(meta)
heatmap_data <- merge(heatmap_data, meta[, c("sample", "Treatment", "Timepoint")], by = "sample")
colnames(heatmap_data)

#heatmap_data_no_Sm <- heatmap_data[, -16]
#colnames(heatmap_data_no_Sm)

dim(heatmap_data)
# Subset the data by Treatment group (Control and Antibiotics)
control_data <- subset(heatmap_data, Treatment == "CONTROL")
antibiotics_data <- subset(heatmap_data, Treatment == "AB")


# Generate heatmap for Control group, sorted by timepoints
generate_heatmap_by_treatment(control_data, "Control")

# Generate heatmap for Antibiotics group, sorted by timepoints
generate_heatmap_by_treatment(antibiotics_data, "Antibiotics")




# Antimicrobial Resistance Analysis

setwd("/Users/jonathan/Documents/Honours/thesis/pipeline/public-data/preterm/amr")
amr_table <- read.csv("summary.tab",header = T,row.names = 1, sep = "\t", check.names=FALSE)
amr_table <- as.data.frame(amr_table)
amr_table$Patient <- sapply(strsplit(sapply(strsplit(rownames(amr_table), "-"), `[`, 1), "/"), `[`, 2)
amr_table
# Load the libraries
library(gt)
library(webshot2)

# Convert your data frame to a gt table
#amr_table_gt <- gt(amr_table)
selected_columns <- amr_table[, c("Patient", "blaSST-1", "aac(6')_Serra", "blaZ", "tet(41)"), drop = FALSE]
amr_table_gt <- gt(selected_columns)
amr_table_gt

# Customize the table (optional)
amr_table_gt <- amr_table_gt %>%
  fmt_number(
    columns = everything(),
    decimals = 2
  ) %>%
  cols_move_to_start(
    columns = "Patient"
  ) %>%
  tab_spanner(
    label = "Genes",
    columns = c(`blaSST-1`, `aac(6')_Serra`, `blaZ`, `tet(41)`)
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
    everything() ~ px(50)
  )

amr_table_gt

#gtsave(amr_table_gt, "amr_table.png", expand = 10)







#### UNUSED CODE


# SPaghetti plots
#library(ggplot2)
#div.taxa<-data.frame(n.otus, n.otus.rar, shann, simp)
#div.taxa
#div.taxa <- rownames_to_column(div.taxa, var = "sample")


#meta
#meta <- rownames_to_column(meta, var = "sample")
#meta
#div.taxa


#spag_data <- merge(div.taxa[, c("sample", "simp")], meta[, c("sample", "Treatment", "Timepoint")], by = 'sample')
#spag_data

#spag_data$sample_id <- sub("-.*", "", spag_data$sample)
#spag_data

#p <- ggplot(data = spag_data, aes(x = Timepoint, y = simp, group = sample_id))
#p + geom_line() + facet_grid(. ~ Treatment)

#p + geom_line() + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
 #                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + facet_grid(. ~ Treatment)
#

# Load RColorBrewer for pastel colors
#library(RColorBrewer)

#colors <- c("lightblue", "lightgreen", "lightpink", "lightyellow", 
 #           "lightcoral", "lightskyblue", "lightgoldenrod", "lightcyan", 
  #          "lightblue", "lightgreen", "lightpink", "lightyellow", 
   #         "lightcoral", "lightskyblue", "lightgoldenrod", "lightcyan")
#
# Create the base plot
#p <- ggplot(data = spag_data, aes(x = Timepoint, y = simp, group = sample_id, color = sample_id))

# Add lines for each sample_id with default R colors
#p + geom_line() + 
 # stat_smooth(aes(group = 1)) + 
  #stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) + 
  #facet_grid(. ~ Treatment) +
  #scale_color_manual(values = colors) +  # Assign default R colors
  #theme_minimal()

#richness

richness.taxa<-rowSums(taxa >0)

meta_ab <- meta[meta$Treatment == "AB",]
meta_ab
#simp_ab <- simp[names(simp) %in% rownames(meta_ab)]
#simp_ab
richness.taxa <- richness.taxa[names(richness.taxa) %in% rownames(meta_ab)]
richness.plot <- data.frame(richness.taxa, meta_ab)
p_value <- kruskal.test(richness.taxa, meta_ab$Timepoint)$p.value


ggplot(richness.plot, aes(x=Timepoint, y=richness.taxa, fill=Timepoint, shape=Timepoint)) + theme_bw() +geom_boxplot() +
  geom_point(position=position_jitterdodge()) + 
  #scale_color_manual(values=c("orange","green","red"))+
  #scale_shape_manual(values=c(21,22)) + 
  #scale_fill_manual(values=c("orange","green","red")) + 
  ylab("Diversity Richness")+ 
  xlab(label = "Treatment") +
  ggtitle(label = "Early Antibiotics") +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.text.x=element_markdown(face="bold",size=10),
        axis.text.y=element_markdown(face="bold",size=14),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
  theme(legend.spacing.y=unit(0.5,"cm"))+
  theme(axis.title.x=element_markdown(face="bold",size=14),
        axis.title.y=element_markdown(face="bold",size=14)) +
  annotate("text", x = 1.5, y = max(richness.plot$richness) + 0.3, 
           label = paste("p-value =", round(p_value, 4)), size = 4, hjust = 0, color="black")


#statistic

wilcox.test(shann~meta$Group)



#Beta diversity..

#------ This is final Filtered Relative Abundance ASV table (normalization)
dim(taxa.filtered)

meta_tp <- meta[meta$Timepoint == "T1",]
meta_tp <- rownames_to_column(meta_tp, var = "sample")
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

