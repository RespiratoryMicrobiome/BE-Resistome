library(phyloseq)
library(vegan)
library(ggplot2)

# Set parameters
wkdir <- file.path("C:/Users/mmaca/Dropbox/Research/Publications/44. CAMEB_2/The_Bronchiectasis_Resistome/Manuscript/CIRCULATED 10.06.2023/Final Submission/Rebuttal/Ivan-Rebuttal/")

ps_AMR <- readRDS(file.path(wkdir, "AMR", "ps_hits.RData"))
ps_Species <- readRDS(file.path(wkdir, "KAIJU", "ps_Species.RData"))
depths <- read.csv(file.path(wkdir, "DEPTHS", "DEPTH-summary2.csv"))
colnames(depths)[1] <- "SampleID"
Master <-read.csv("../Data/R_input_files//Clinical_AMR_Microbiome_R2.csv") %>%
  as_tibble()
MasterVIZ = Master
MasterVIZ$select <- ifelse(MasterVIZ$SC_AMR_alt==0, "null", "Bronchiectasis")
MasterVIZ$select <- ifelse(is.na(MasterVIZ$select), "Non-diseased", MasterVIZ$select)
MasterVIZ$SC_AMR_alt <- ifelse(is.na(MasterVIZ$SC_AMR_alt), "Non-diseased", MasterVIZ$SC_AMR_alt)
AMRDiversityViz<-subset(MasterVIZ, select != "null")

N_Species <- apply(otu_table(ps_Species), 1, function(x) sum(x!=0))
N_AMRgenes <- apply(otu_table(ps_AMR), 1, function(x) sum(x!=0))
N_PAreads <- c(otu_table(ps_Species)[, "Pseudomonas aeruginosa"])
SpeciesDiversity <- diversity(otu_table(ps_Species), index= "shannon")
PA <- sapply(otu_table(ps_Species)[, "Pseudomonas aeruginosa"], function(x) if(x>50) "Pos" else "Neg")
KP <- sapply(otu_table(ps_Species)[, "Klebsiella pneumoniae"], function(x) if(x>50) "Pos" else "Neg")
PAKP <- paste0(PA, ".", KP)
dfSpecies <- data.frame(SampleID= sample_names(ps_Species), N_Species= N_Species,
  PA= PA, KP= KP, PAKP= PAKP, N_PAreads= N_PAreads, SpeciesDiversity= SpeciesDiversity)
dfSpecies$PAKP <- 
  gsub("Neg.Neg", "PA-KP-", gsub("Neg.Pos", "PA-KP+", gsub("Pos.Neg", "PA+KP-", gsub("Pos.Pos", "PA+KP+", dfSpecies$PAKP))))
AMRDiversity <- diversity(otu_table(ps_AMR), index= "shannon")
dfAMR <- data.frame(SampleID= sample_names(ps_AMR), N_AMRgenes= N_AMRgenes, AMRDiversity= AMRDiversity)

df <- depths
df <- merge(x= df, y= dfSpecies, by= "SampleID", all.x= TRUE)
df <- merge(x= df, y= dfAMR, by= "SampleID", all.x= TRUE)
rownames(df) <- df$SampleID
df$NonhostDepth <- sapply(df$N_nonhuman, function(x) if (x>100000) "Deep" else "Shallow")
table(df$NonhostDepth)
df$NonhostDepth <- factor(df$NonhostDepth, levels = c("Shallow", "Deep"))


#CHECK
#all(sample_names(ps_AMR)==sample_names(ps_Species))
#all(rownames(df)==sample_names(ps_Species))
sample_data(ps_AMR) <- df
sample_data(ps_Species) <- df

df <- merge(df, MasterVIZ[, c("SampleID", "SC_AMR_alt")], by = "SampleID", all.x = TRUE)

boxplot(N_nonhuman ~ SC_AMR_alt, data= df[which(df$SC_AMR_alt != 0 ),])
wilcox.test(N_nonhuman ~ SC_AMR_alt, data= df[which(df$SC_AMR_alt != 0 ),])

boxplot(Percent_nonhuman ~ SC_AMR_alt, data= df[which(df$SC_AMR_alt != 0 ),])
wilcox.test(Percent_nonhuman ~ SC_AMR_alt, data= df[which(df$SC_AMR_alt != 0 ),])

boxplot(Percent_nonhuman ~ SC_AMR_alt, data= df[which(df$SC_AMR_alt != 0 ),])
wilcox.test(Percent_nonhuman ~ SC_AMR_alt, data= df[which(df$SC_AMR_alt != 0 ),])

boxplot(N_Species ~ NonhostDepth, data= df)
wilcox.test(N_Species ~ NonhostDepth, data= df)

chisq.test(df[which(df$SC_AMR_alt != 0 ),]$NonhostDepth, df[which(df$SC_AMR_alt != 0 ),]$SC_AMR_alt)

table(df$NonhostDepth, df$SC_AMR_alt)

boxplot(SpeciesDiversity ~ NonhostDepth, data= df)
wilcox.test(SpeciesDiversity ~ NonhostDepth, data= df)

plot(df$N_nonhuman, df$SpeciesDiversity, xlab= "No. of non-host reads", ylab= "Species diversity")
cor.test(df$N_nonhuman, df$SpeciesDiversity)

plot(df0$Percent_nonhuman, df0$AMRDiversity, xlab= "No. of non-host reads", ylab= "AMR diversity")
cor.test(df0$Percent_nonhuman, df0$AMRDiversity)

df0 <- df[df$AMRDiversity!=0, ]

# Perform the correlation test and store results
cor_test_result.a <- cor.test(log10(df0$N_nonhuman), df0$AMRDiversity)

# Extract the correlation coefficient and p-value
cor_coefficient.a <- cor_test_result.a$estimate
p_value.a <- cor_test_result.a$p.value

Pt2<-ggplot(df0, aes(x = N_nonhuman, y = AMRDiversity)) +
  geom_point(aes(color = NonhostDepth), alpha = 0.6) +  # Color grouping only for points
  geom_smooth(method = "loess", se = TRUE, color = "black") +  # Overall regression line
  labs(title = "Non-Human reads vs ARG diversity",
       x = "Non-Human reads",
       y = "ARG diversity (SDI)",
       color = "Sequencing depth") +  # Rename the color legend
  scale_x_log10()+
  theme_minimal(base_size = 14) +
  annotate("text", x = quantile(df0$N_nonhuman, 0.85), 
           y = quantile(df0$AMRDiversity, 0), 
           label = sprintf("r = %.2f, p = %.3f", cor_coefficient.a, p_value.a), 
           size = 5, hjust = 0)

# Perform the correlation test and store results
cor_test_result.s <- cor.test(df0$N_nonhuman, df0$SpeciesDiversity)

# Extract the correlation coefficient and p-value
cor_coefficient.s <- cor_test_result.s$estimate
p_value.s <- cor_test_result.s$p.value

Pt1<-ggplot(df0, aes(x = N_nonhuman, y = SpeciesDiversity)) +
  geom_point(aes(color = NonhostDepth), alpha = 0.6) +  # Color grouping only for points
  geom_smooth(method = "loess", se = TRUE, color = "black") +  # Overall regression line
  labs(title = "Non-Human reads vs Species diversity",
       x = "Non-Human reads",
       y = "Species diversity (SDI)",
       color = "Sequencing depth") +  # Rename the color legend
  theme_minimal(base_size = 14) +
  annotate("text", x = quantile(df0$N_nonhuman, 0.85), 
           y = quantile(df0$SpeciesDiversity, 0), 
           label = sprintf("r = %.2f, p = %.3f", cor_coefficient.s, p_value.s), 
           size = 5, hjust = 0) +
  scale_x_log10() +
  guides(color = guide_legend(reverse = TRUE))

df0 <- df[df$SC_AMR_alt!=0, ]

# Perform Wilcoxon test
wilcox_result <- wilcox.test(N_nonhuman ~ SC_AMR_alt, data = df0)

# Create the plot
Pt3 <- ggplot(data = df0, aes(x = SC_AMR_alt, y = N_nonhuman, group = SC_AMR_alt)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(fill = SC_AMR_alt)) +
  geom_jitter(alpha = 0.6, width = 0.2, aes(color = SC_AMR_alt)) +
  scale_fill_manual(values = c("#1800F5", "#932DE7")) +
  scale_color_manual(values = c("#1800F5", "#932DE7")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("RT1", "RT2")) +
  labs(title = "Non-Human reads RT1 vs RT2",
       x = "Resistotype",
       y = "log10(non-human reads)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  annotate("text", x = 1.5, y = max(df0$N_nonhuman),  # Adjust x and y for positioning
           label = sprintf("Wilcoxon p-value: %.3f", wilcox_result$p.value),
           size = 4, vjust = 1)  # Adjust text size and vertical position

wilcox_result.pc <- wilcox.test(Percent_nonhuman ~ SC_AMR_alt, data = df0)

# Create the plot
Pt4 <- ggplot(data = df0, aes(x = SC_AMR_alt, y = Percent_nonhuman, group = SC_AMR_alt)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(fill = SC_AMR_alt)) +
  geom_jitter(alpha = 0.6, width = 0.2, aes(color = SC_AMR_alt)) +
  scale_fill_manual(values = c("#1800F5", "#932DE7")) +
  scale_color_manual(values = c("#1800F5", "#932DE7")) +
  #scale_y_log10() +
  scale_x_discrete(labels = c("RT1", "RT2")) +
  labs(title = "Non-Human reads RT1 vs RT2",
       x = "Resistotype",
       y = "Non-human reads(%)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  annotate("text", x = 1.5, y = max(df0$Percent_nonhuman),  # Adjust x and y for positioning
           label = sprintf("Wilcoxon p-value: %.3f", wilcox_result.pc$p.value),
           size = 4, vjust = 1)  # Adjust text size and vertical position


# Assuming Pt1, Pt2, and Pt3 are your ggplot objects
FigureSX<-ggarrange(Pt1, Pt2, Pt3, Pt4,
          ncol = 2,   # Arrange in 3 columns
          nrow = 2,   # Arrange in 1 row
          common.legend = TRUE, # Use a common legend
          legend = "bottom")    # Place the legend at the bottom

#ggsave("C:/Users/mmaca/Code/git/BE-Resistome/Data/R_output_files/Pt3_plot.pdf", plot = FigureSX, width = 10, height = 8)


# Set parameters
wkdir <- file.path("C:/Users/mmaca/Dropbox/Research/Publications/44. CAMEB_2/The_Bronchiectasis_Resistome/Manuscript/CIRCULATED 10.06.2023/Final Submission/Rebuttal/Ivan-Rebuttal/")

ps_AMR <- readRDS(file.path(wkdir, "AMR", "ps_hits.RData"))
ps_Species <- readRDS(file.path(wkdir, "KAIJU", "ps_Species.RData"))
depths <- read.csv(file.path(wkdir, "DEPTHS", "DEPTH-summary2.csv"))
colnames(depths)[1] <- "SampleID"

N_Species <- apply(otu_table(ps_Species), 1, function(x) sum(x!=0))
N_AMRgenes <- apply(otu_table(ps_AMR), 1, function(x) sum(x!=0))

ps_Species.prop <- transform_sample_counts(ps_Species, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
PA <- sapply(otu_table(ps_Species.prop)[, "Pseudomonas aeruginosa"], function(x) if(x>50) ">50" else "<=50")
dfSpecies <- data.frame(SampleID= sample_names(ps_Species), N_Species= N_Species, PA= PA)
dfAMR <- data.frame(SampleID= sample_names(ps_AMR), N_AMRgenes= N_AMRgenes)

df <- depths
df <- merge(x= df, y= dfSpecies, by= "SampleID", all.x= TRUE)
df <- merge(x= df, y= dfAMR, by= "SampleID", all.x= TRUE)

boxplot(N_Species ~ PA, data= df)
wilcox.test(N_Species ~ PA, data= df)
tapply(df$N_Species, df$PA, summary)
tapply(df$N_Species, df$PA, IQR)

boxplot(N_nonhuman ~ PA, data= df)
wilcox.test(N_nonhuman ~ PA, data= df)
tapply(df$N_nonhuman, df$PA, summary)
tapply(df$N_nonhuman, df$PA, IQR)