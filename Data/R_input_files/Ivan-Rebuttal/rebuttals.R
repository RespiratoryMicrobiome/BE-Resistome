library(phyloseq)
library(vegan)

# Set parameters
wkdir <- file.path("C:/Users/fivan/Downloads/CAMEB2 Resistome/BE-Resistome-main")

ps_AMR <- readRDS(file.path(wkdir, "AMR", "ps_hits.RData"))
ps_Species <- readRDS(file.path(wkdir, "KAIJU", "ps_Species.RData"))
depths <- read.csv(file.path(wkdir, "DEPTHS", "DEPTH-summary2.csv"))
colnames(depths)[1] <- "SampleID"

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

#CHECK
#all(sample_names(ps_AMR)==sample_names(ps_Species))
#all(rownames(df)==sample_names(ps_Species))
sample_data(ps_AMR) <- df
sample_data(ps_Species) <- df

boxplot(N_AMRgenes ~ NonhostDepth, data= df)
wilcox.test(N_AMRgenes ~ NonhostDepth, data= df)

boxplot(N_Species ~ NonhostDepth, data= df)
wilcox.test(N_Species ~ NonhostDepth, data= df)

boxplot(AMRDiversity ~ NonhostDepth, data= df)
wilcox.test(AMRDiversity ~ NonhostDepth, data= df)

boxplot(SpeciesDiversity ~ NonhostDepth, data= df)
wilcox.test(SpeciesDiversity ~ NonhostDepth, data= df)

adonis2(as.data.frame(otu_table(ps_Species)) ~ NonhostDepth, data= as.data.frame(as.matrix(sample_data(ps_Species))))
adonis2(as.data.frame(otu_table(ps_AMR)) ~ NonhostDepth, data= as.data.frame(as.matrix(sample_data(ps_AMR))))

library(lattice)
with(df, xyplot(N_AMRgenes ~ N_nonhuman, xlab= "No. of non-host reads", ylab= "No. of resistance genes",
  group= PA, auto.key= list(space = "right")))
with(df, xyplot(N_Species ~ N_nonhuman, xlab= "No. of non-host reads", ylab= "No. of species",
  group= PA, auto.key= list(space = "right")))

boxplot(N_Species ~ PA, data= df)
wilcox.test(N_Species ~ PA, data= df)

plot(df$N_PAreads, df$Percent_AMR, xlab= "No. of PA reads", ylab= "Proportion of AMR hits")
cor.test(df$N_PAreads, df$Percent_AMR)

plot(df$N_nonhuman, df$SpeciesDiversity, xlab= "No. of non-host reads", ylab= "Species diversity")
cor.test(df$N_nonhuman, df$SpeciesDiversity)

df0 <- df[df$AMRDiversity!=0, ]
plot(df0$N_nonhuman, df0$AMRDiversity, xlab= "No. of non-host reads", ylab= "AMR diversity")
cor.test(df0$N_nonhuman, df0$AMRDiversity)





library(phyloseq)
library(vegan)

# Set parameters
wkdir <- file.path("C:/Users/fivan/Downloads/CAMEB2 Resistome/BE-Resistome-main")

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











