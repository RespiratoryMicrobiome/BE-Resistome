#SETUP####
#Set the working directory
#setwd("/Users/michealmacaogain/Dropbox/Research/Publications/44. CAMEB_2/Take_over/MASTER DATA/")

#setwd("/home/mmacaogain/Dropbox/Research/Publications/44. CAMEB_2/Take_over/MASTER DATA/")

setwd("C:/Users/mmaca/Dropbox/Research/Publications/44. CAMEB_2/Take_over/MASTER DATA/git_repo/")

#Load packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, ggplot2, tidyverse, tidyr, lubridate, data.table, tsibble, wesanderson, reticulate, SNFtool, funrar, vegan, dunn.test, ggpubr, Hmisc, RColorBrewer, phyloseq, dplyr, reshape2, forcats, colorspace, pheatmap)

#reticulate env set up for python. 
reticulate::use_python("C:/Users/mmaca/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe", required = TRUE)
source("R_input_files/function_snf.R")
source_python("R_input_files/sil.py")

# Assuming your dissimilarity matrix is ab_dsim
# and you want to check number of clusters from 2 to 10

install.packages('fpc')
library('fpc')

# Initialize a vector to store average silhouette widths
avg_sil_widths <- c()

for (k in 2:10) {
  
  # Perform clustering with k clusters
  # Replace `spectralClustering` with the appropriate clustering function from the `SNFtool` package
  labels <- spectralClustering(AB, k)
  
  # Compute cluster statistics
  clust_stats <- cluster.stats(ab_dsim, labels)
  
  # Store average silhouette width
  avg_sil_widths <- c(avg_sil_widths, clust_stats$avg.silwidth)
}

# Print average silhouette widths
print(avg_sil_widths)


#DATA####
##Master data####
Master <-read.csv("R_input_files/Clinical_AMR_Microbiome.csv") %>%
  as_tibble()
Master$FEVfactor<-cut(Master$FEV1, breaks=c(0, 30, 50, 70, Inf))

##Longitudinal AMR data ####
MasterLT <-read.csv("R_input_files/LT_master_combined_8.0.csv")

###wrangle AMR data ####
AMRFam <- Master %>% #clinical variables + amr families
  as_tibble() %>%
  select(-29:-42,-64:-356)
AMRFam$FEVfactor<-cut(AMRFam$FEV1,  breaks=c(0, 30, 50, 70, Inf))
#set levels
AMRFam$ExacerbatorState <- factor(AMRFam$ExacerbatorState, levels=c("NonEx", "Exacerbator", "FreqEx"))
AMRFam$Country <- factor(AMRFam$Country, levels=c("SG", "KL", "DD", "MI"))
AMRFam$Aetiology_short <- factor(AMRFam$Aetiology_short, levels=c("idiopathic", "postInfect", "postTB", "other"))
AMRFam$SampleID <- factor(AMRFam$SampleID, levels = AMRFam$SampleID[order(AMRFam$SC_AMR_alt)])
AMRFam$FEVfactor<-fct_rev(AMRFam$FEVfactor)
AMRFam <- AMRFam %>% 
  gather(Resistome, RPKM, Acridine.dye,Aminocoumarin.antibiotic,Aminoglycoside,Antibacterial.free.fatty.acids,Beta.lactam,Bicyclomycin,Diaminopyrimidine,Fluoroquinolone,Fosfomycin,Fusidic.acid,MLS,Multidrug,Mupirocin,Nitroimidazole.antibiotic,Nucleoside.antibiotic,Peptide.antibiotic,Phenicol,Rifampicin,Sulfonamide.antibiotic,Tetracycline,Triclosan, -SampleID, -Country,  -Continent, -Matching, -Paired, -Trio, -Age,-Sex..Male.0..Female.1.,-Exacerbations,-ExacerbatorState,  -FEV1,   -BSI, -ICS.use,   -BMI, -Aetiology, -Aetiology_short,-MMRC.score,-SC_AMR_alt, -FEVfactor)

AMRFam$CTRL<-ifelse(is.na(AMRFam$Age), "CTRL", "PATIENT")

###wrangle Bphage data####
BphageFam <- Master %>%
  as_tibble() %>%
  select(-43:-355)
#filter(ExacerbatorState != "NA") #%>%
#filter(Continent == "Asia") %>%
#filter(Matching == "Matched" )
#generate a % value to rank by myoviridae / Phycodnaviridae
IridoviridaePC<-(BphageFam$Iridoviridae/(rowSums(BphageFam[, c(28:41)]))*100)
BphageFam$Country <- factor(BphageFam$Country, levels=c("SG", "KL", "DD", "MI"))
BphageFam$SampleID <- factor(BphageFam$SampleID, levels = rev(BphageFam$SampleID[order(IridoviridaePC)]))

#gather on data
BphageFam <- BphageFam %>% 
  gather(Virome, RPKM, Siphoviridae, Unassigned, Iridoviridae, Myoviridae, Phycodnaviridae, Polydnaviridae, Picornaviridae, Podoviridae, Poxviridae, Nyamiviridae, Mimiviridae, Herpesviridae, Inoviridae, Alloherpesviridae, -SampleID, -Country,  -Continent, -Matching, -Paired, -Trio, -Age,-Sex..Male.0..Female.1.,-Exacerbations,-ExacerbatorState,  -FEV1,   -BSI, -ICS.use,   -BMI, -Aetiology, -Aetiology_short,-MMRC.score,-SC_AMR_alt, -FEVfactor)
#post gather level setting (required??)
BphageFam$Virome <- factor(BphageFam$Virome, levels = c("Iridoviridae", "Siphoviridae","Myoviridae", "Phycodnaviridae","Polydnaviridae","Picornaviridae", "Poxviridae","Podoviridae","Nyamiviridae","Mimiviridae","Herpesviridae","Inoviridae","Alloherpesviridae", "Unassigned"))

#Spectral clustering####
##AMR data lables####
ab_data=make_relative(as.matrix(read.csv("../Spectral_clustering/AMR.csv", row.names = 1)))*100 #load functional data
ab_data[is.nan(ab_data)] <- 0
ab_data<-as.data.frame(ab_data)
#filter based on prevalence
z=colSums(ab_data>0.1) #filter to reduced number of taxa 
sel_col=row.names(as.data.frame(z[z>=(0.01*(nrow(ab_data)))])) #In 1% patients prevalent
ab_data<-ab_data[sel_col]
remove(sel_col,z)
ab_data<-ab_data[rowSums(ab_data[, -1])>0, ] #drop no_res samples
ab_data<-ab_data[row.names(ab_data) != "TBS672", , drop = FALSE] #TBS672 contains a single PatA gene creating clustering artifact

ab_data<-as.data.frame(make_relative(as.matrix(ab_data)))
#write.csv(ab_data, "R_output_files/AMRFiltered.csv")

#Spectral clustering
#create vegdist similarity matrix
ab_dsim=vegdist(ab_data,method='bray',diag=TRUE,upper=TRUE) 
ab_dsim[is.nan(ab_dsim)]<-0 
AB=(as.matrix(ab_dsim)-1)*-1
heatmap(as.matrix(ab_dsim))
#tune for K
sil_values=c()
for (i in 2:20){
  labels = spectralClustering(AB, i)
  sil_values = c(sil_values, silhouette_score(AB, labels))
}
tuned_k<-which.max(sil_values)
#assess clusters and grouping
paste(tuned_k,sil_values[tuned_k],sep = " ")
labels=spectralClustering(AB,tuned_k)
labels <- max(labels)+1 - labels #aesthetic change: most prevalent label now assigned value = 1. 
table(labels)
lab=as.data.frame(labels,row.names = row.names(ab_data))

#Figure 1 ####
##Stacked Barplot AMR rel abudance ####
HvsBE<-ggplot(data=AMRFam,aes(x=CTRL, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('Non-diseased','Bronchiectasis'))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(ncol=1), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

## PCA plot Non-diseased vs Bronchiectasis####
AMR_diversity <- Master %>%
  as_tibble() %>%
  #select(1:1,395:645) #for genes
  select(1:1,43:63) #for amr drug class
NAMES_list <- AMR_diversity$SampleID
main_data <- AMR_diversity[AMR_diversity$SampleID %in% NAMES_list, ]
AMR_diversity<-as.matrix(AMR_diversity)
rownames(AMR_diversity) <- AMR_diversity[,1]
AMR_diversity = as.data.frame(subset(AMR_diversity, select = -c(SampleID) ))
AMR_diversity[] <- lapply(AMR_diversity, as.numeric)
AMR_diversity<-AMR_diversity[row.names(AMR_diversity) != "TBS672", , drop = FALSE]

isZero <- base::rowSums(AMR_diversity) == 0
#sum(isZero)#NO amr detected in 37 samples
AMR_diversity<-AMR_diversity[!isZero,]

vegdist(AMR_diversity, "bray")-> Mbiome_PCoA
as.matrix(Mbiome_PCoA)->Mbiome_PCoA
BrayCurtMbiome=cmdscale(Mbiome_PCoA)
#ordiplot (BrayCurtMbiome, display = 'species', type = 'text')
BCords<-scores(BrayCurtMbiome)
BCords<-(as.data.frame(t(BCords)))
BCords<-as.data.frame(t(BCords))

MasterVIZ = Master
MasterVIZ$select <- ifelse(MasterVIZ$SC_AMR_alt==0, "null", "Bronchiectasis")
MasterVIZ$select <- ifelse(is.na(MasterVIZ$select), "Non-diseased", MasterVIZ$select)
MasterVIZ$SC_AMR_alt <- ifelse(is.na(MasterVIZ$SC_AMR_alt), "Non-diseased", MasterVIZ$SC_AMR_alt)
AMRDiversityViz<-subset(MasterVIZ, select != "null")

AMRDiversityViz<-AMRDiversityViz[AMRDiversityViz$SampleID != "TBS153", , drop = FALSE] #remove for gene level analysis

AMRDiversityViz$Dim1<-BCords$Dim1
AMRDiversityViz$Dim2<-BCords$Dim2

#checking PC %s
#rda(X = AMR_diversity, scale = TRUE)
#assess explained variation
checkEig<-capscale(AMR_diversity ~1)
Eig <-eigenvals(checkEig)
Eig / sum(Eig)

#AMR PCOA of Resistotypes BY SC_RESISTOTYPE   
gg <- data.frame(cluster=factor(AMRDiversityViz$select), x=AMRDiversityViz$Dim1, y=AMRDiversityViz$Dim2, grp=AMRDiversityViz$select)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
BC<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_shape_discrete(labels = c("Healthy", "Bronchiectasis"))+
  scale_colour_manual(values = c("#F8766D", "#619CFF"), labels = c("Bronchiectasis","Healthy"))+
  labs(colour="",  
       x = "PC 1 (77.6%)", y = "PC 2 (6.9%)")+
  theme(legend.position="none",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_x_reverse()+
  #scale_y_reverse()+ #add for gene level analysis
  guides(colour = guide_legend(reverse = T))

#Are patients distinct from healthy controls?
adonis2(AMR_diversity ~ select, data=AMRDiversityViz, method="bray", permutations=9999)

##Longitudinal data####
###Data wrangle####
AMRLT <- MasterLT %>%
  as_tibble() %>%
  select(-14:-230) #Version 7.0 has now lost the bacteriophage analysis columns consider dropping these. Have created version 8.0 dropping bacteriophage data
AMR_cols<-colnames(AMRLT[14:34])
AMRLT <- AMRLT %>% 
  gather(AMR, RPKM, AMR_cols, -SampleSeqNo, -SputumSampleNo,  -TypeSamples, -TypeSamplesA,-TypeSamplesB,-Exacerbations,	-FEV1, -BSI,	-Severity, -WksToNxtEx,	-TmToNxtEx,	-Antibiotic,	-Antibiotic_class)
AMRLT$TmToNxtEx <- factor(AMRLT$TmToNxtEx , levels = c("MoreThan12w","LessThan12w"))
AMRLT$Exacerbations <- factor(AMRLT$Exacerbations , levels = c("NFE","FE"))
relapse.labs <- c(
  `LessThan12w` = "<12 w",
  `MoreThan12w` = ">12 w")
AMRLT$FEV170<-ifelse(AMRLT$FEV1 >70, ">70", "<70")
AMRLTctrols<-subset(AMRLT, is.na(TypeSamplesB))

###Barplot stable####
STBL<-ggplot(data=AMRLTctrols,aes(x=SampleSeqNo, y=RPKM, fill=AMR))+
  geom_bar(aes(), stat="identity", position = 'fill') +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_x_discrete(labels = c('BL','FUP'))+
  scale_y_continuous(labels = scales::percent)+
  theme(#legend.position="none",
    #axis.text=element_blank(),
    axis.title=element_text(size=14),
    #axis.text.x = element_text(angle = 90),
    panel.background = element_rect(fill = NA),
    axis.line = element_line(size = 0.5, colour = "black"))+
  #legend.title = element_blank(),
  #legend.text = element_text(face = "italic"))
  #guides(fill=guide_legend(ncol=1), size = 0.1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRLTctrols$SputumSampleNo, scales="free_x")+
  theme(strip.background = element_rect(
    color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    #legend.text=element_text(size=8)
    legend.position="none"
  )

#drop the controls stable
AMRLT<-subset(AMRLT, is.na(Severity) != TRUE)

###Barplots Exacerbation####
C<-ggplot(data=AMRLT[which(AMRLT$TypeSamplesA !="NA"),],aes(x=TypeSamplesB, y=RPKM, fill=AMR))+
  geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('B','E', 'P'))+
  theme(#legend.position="none",
    #axis.text=element_blank(),
    axis.title=element_text(size=14),
    #axis.text.x = element_text(angle = 90),
    panel.background = element_rect(fill = NA),
    axis.line = element_line(size = 0.5, colour = "black"))+
  #legend.title = element_blank(),
  #legend.text = element_text(face = "italic"))
  #guides(fill=guide_legend(ncol=1), size = 0.1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_grid(~AMRLT$TmToNxtEx, scales="free_x", labeller = as_labeller(relapse.labs))+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    #legend.text=element_text(size=8)
    legend.position="none"
  )

E<-ggplot(data=AMRLT[which(AMRLT$TypeSamplesA !="NA"),],aes(x=TypeSamplesB, y=RPKM, fill=AMR))+
  geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('B','E', 'P'))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"))+
  #legend.title = element_blank(),
  #legend.text = element_text(face = "italic"))
  #guides(fill=guide_legend(ncol=1), size = 0.1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRLT[which(AMRLT$TypeSamplesA !="NA"),]$TmToNxtEx, scales="free_x", labeller = as_labeller(relapse.labs))+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    #legend.text=element_text(size=8)
    legend.position="none"
  )

FEV<-ggplot(data=AMRLT[which(AMRLT$TypeSamplesA !="NA"),],aes(x=TypeSamplesB, y=RPKM, fill=AMR))+
  geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('B','E', 'P'))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"))+
  #legend.title = element_blank(),
  #legend.text = element_text(face = "italic"))
  #guides(fill=guide_legend(ncol=1), size = 0.1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_wrap(~AMRLT[which(AMRLT$TypeSamplesA !="NA"),]$FEV170, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    #legend.text=element_text(size=8)
    legend.position="none"
  )

### PCA ####
AMRLT_diversity <- MasterLT[which(MasterLT$TypeSamplesA !="NA"),] %>%
  as_tibble() %>%
  select(-14:-230)
#AMRLT_diversity <- select(MasterLT[which(MasterLT$TypeSamplesA !="NA"),], -2:-79, -231:-3621) #ugly 'which subsetting to drop controls #gene
AMRLT_diversity <- select(MasterLT[which(MasterLT$TypeSamplesA !="NA"),], -2:-230) #ugly 'which subsetting to drop controls #class

NAMES_list <- head(MasterLT$SampleSeqNo, -6) #head is just to drop controls again n=6
main_dataLT <- AMRLT_diversity[AMRLT_diversity$SampleSeqNo %in% NAMES_list, ]
AMRLT_diversity<-as.matrix(AMRLT_diversity)
rownames(AMRLT_diversity) <- AMRLT_diversity[,1]
AMRLT_diversity = as.data.frame(subset(AMRLT_diversity, select = -c(SampleSeqNo) ))
AMRLT_diversity[] <- lapply(AMRLT_diversity, as.numeric)
isZero <- base::rowSums(AMRLT_diversity) == 0
sum(isZero)
AMRLT_diversity<-AMRLT_diversity[!isZero,]
vegdist(AMRLT_diversity, "bray")-> Mbiome_PCoA
as.matrix(Mbiome_PCoA)->Mbiome_PCoA
BrayCurtMbiome=cmdscale(Mbiome_PCoA)
#ordiplot (BrayCurtMbiome, display = 'species', type = 'text')
BCords<-scores(BrayCurtMbiome)
BCords<-(as.data.frame(t(BCords)))
BCords<-as.data.frame(t(BCords))

LTDiversityViz<-MasterLT[which(MasterLT$TypeSamplesA !="NA"),] #drop controls -which subsetting
#LTDiversityViz$SampleSeqNo %in% row.names(BCords) 
LTDiversityViz<-LTDiversityViz[ LTDiversityViz$SampleSeqNo %in% row.names(BCords) , ]

LTDiversityViz$Dim1<-BCords$Dim1
LTDiversityViz$Dim2<-BCords$Dim2

LTDiversityViz$FEV170<-ifelse(LTDiversityViz$FEV1 >70, ">70", "<70")
LTDiversityViz$FEV170<- factor(LTDiversityViz$FEV170 , levels = c(">70","<70"))

#AMR PCOA of Resistotypes BY sample type   
gg <- data.frame(cluster=factor(LTDiversityViz$TypeSamplesB), x=LTDiversityViz$Dim1, y=LTDiversityViz$Dim2, grp=LTDiversityViz$TypeSamplesB, shape=LTDiversityViz$TypeSamplesB)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
D<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_shape_discrete(labels = c("B", "E", "P"))+
  scale_colour_manual(values = c("#619CFF", "#F8766D", "#00BA38"), labels = c("B", "E", "P"))+
  labs(colour="",  
       x = "PC 1 (83.3%)", y = "PC 2 (12.1%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+  scale_x_reverse()
#+ggtitle("Timepoint")

#AMR PCOA of Resistotypes BY SC_RESISTOTYPE   
gg <- data.frame(cluster=factor(LTDiversityViz$TmToNxtEx), x=LTDiversityViz$Dim1, y=LTDiversityViz$Dim2, grp=LTDiversityViz$TmToNxtEx, shape=LTDiversityViz$TypeSamplesB)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
F<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_shape_discrete(labels = c("B", "E", "P"))+
  scale_colour_manual(values = c("#F8766D", "#619CFF"), labels = c("<12 w", ">12 w"))+
  labs(colour="",  
       x = "PC 1 (83.3%)", y = "PC 2 (12.1%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_x_reverse()

#+ggtitle("Time to next exacerbation")

#AMR PCOA of Resistotypes BY FEV   
gg <- data.frame(cluster=factor(LTDiversityViz$FEV170), x=LTDiversityViz$Dim1, y=LTDiversityViz$Dim2, grp=LTDiversityViz$FEV170, shape=LTDiversityViz$TypeSamplesB)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
G<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_shape_discrete(labels = c("B", "E", "P"))+
  scale_colour_manual(values = c("#619CFF","#F8766D"))+
  labs(colour="",  
       x = "PC 1 (83.3%)", y = "PC 2 (12.1%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),)+
  scale_x_reverse()

#checking PC %s
#rda(X = AMRLT_diversity, scale = TRUE)
checkEig<-capscale(AMRLT_diversity ~1)
Eig <-eigenvals(checkEig)
Eig / sum(Eig)

#does timepoint matter?
adonis2(AMRLT_diversity~TypeSamplesB , data = LTDiversityViz, method = "bray",permutations=9999)
#does time to next exacerbation matter?
adonis2(AMRLT_diversity~TmToNxtEx , data = LTDiversityViz, method = "bray",permutations=9999)

##Combine and print panels for Figure 1####
Figure_1top<-ggarrange(HvsBE,BC,STBL,
                   font.label = list(size = 5),
                   common.legend = FALSE, nrow = 1, ncol = 3) #this one

Figure_1bot<-ggarrange(C,NULL,E,D,NULL,F, font.label = list(size = 5),
                       common.legend = FALSE, widths = c(1, 0.1,1,1,0.1,1))

Figure_1<- ggarrange(Figure_1top, Figure_1bot, font.label = list(size = 5),
                           common.legend = FALSE,heights = c(0.35, 0.65), nrow =2)

pdf(file = "R_output_files/Figures/Figure_1.pdf",   # The directory you want to save the file in
   width = 10, # The width of the plot in inches
  height = 12)
Figure_1
dev.off()

##supplementary ####
S2A<-ggplot(data=AMRLT[which(is.na(AMRLT$TypeSamplesB)==FALSE),],aes(x=TypeSamplesB, y=RPKM, fill=AMR))+ #drop controls...
  geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('B','E', 'P'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"))+
  #legend.title = element_blank(),
  #legend.text = element_text(face = "italic"))
  guides(fill=guide_legend(nrow=3), size = 0.1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRLT[which(is.na(AMRLT$TypeSamplesB)==FALSE),]$Exacerbations, scales="free_x", labeller = as_labeller(relapse.labs))+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    #legend.text=element_text(size=8)
    #legend.position="bottom"
  )
#S2A
S2B<-ggplot(data=AMRLT[which(is.na(AMRLT$TypeSamplesB)==FALSE),],aes(x=TypeSamplesB, y=RPKM, fill=AMR))+
  geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('B','E', 'P'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"))+
  #legend.title = element_blank(),
  #legend.text = element_text(face = "italic"))
  guides(fill=guide_legend(nrow=3), size = 0.1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRLT[which(is.na(AMRLT$TypeSamplesB)==FALSE),]$Severity, scales="free_x", labeller = as_labeller(relapse.labs))+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    #legend.text=element_text(size=8)
    #legend.position="none"
  )
#S2B

#AMR PCOA of Resistotypes LT Exacerbation  
gg <- data.frame(cluster=factor(LTDiversityViz$Exacerbations), x=LTDiversityViz$Dim1, y=LTDiversityViz$Dim2, grp=LTDiversityViz$Exacerbations, shape=LTDiversityViz$TypeSamplesB)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
S2C<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_shape_discrete(labels = c("B", "E", "P"))+
  scale_colour_manual(values = c("#F8766D", "#619CFF"), labels = c("F.Ex", "Non-F.Ex"))+
  labs(colour="",  
       x = "PC 1 (83.3%)", y = "PC 2 (12.1%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+  scale_x_reverse()

#AMR PCOA of Resistotypes LT Severity  
gg <- data.frame(cluster=factor(LTDiversityViz$Severity), x=LTDiversityViz$Dim1, y=LTDiversityViz$Dim2, grp=LTDiversityViz$Severity, shape=LTDiversityViz$TypeSamplesB)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
S2D<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_shape_discrete(labels = c("B", "E", "P"))+
  scale_colour_manual(values = c("#619CFF", "#00BA38","#F8766D"), labels = c("Mild", "Moderate", "Severe"))+
  labs(colour="",  
       x = "PC 1 (83.3%)", y = "PC 2 (12.1%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+  scale_x_reverse()
Figure_S2_temp1<-ggarrange(S2A, S2B,
                     font.label = list(size = 5),
                     common.legend = TRUE, legend = "bottom")
Figure_S2_temp2<-ggarrange(S2C, S2D,
                    font.label = list(size = 5),
                    common.legend = FALSE, nrow = 2, ncol = 2) 
Figure_S2<-ggarrange(S2A, S2B,S2C, S2D,
                     font.label = list(size = 5),common.legend = TRUE )

#this one
pdf(file = "R_output_files/Figures/Figure_S2.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12)
Figure_S2
dev.off()

#Blank Analysis####
BlankData <-read.csv("./R_input_files/Blanks/blank_analysis.csv") %>%
  as_tibble() #%>%   

TaxaBlank <- BlankData  %>%
  select(1:42)
amrBlank <- BlankData  %>%
  select(1:2,43:66)

TaxaBlank<-melt(TaxaBlank, id.vars = c("Sample_ID", "Type"))
TaxaBlank$Type<-factor(TaxaBlank$Type, levels= c("Sample","Blank","Blank-seq"))

amrBlank<-melt(amrBlank, id.vars = c("Sample_ID", "Type"))

n <- 41
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_spec<-replace(col_vector, 41, "grey") #

Blank_taxa<-ggplot(data=TaxaBlank,aes(x=Type, y=value, fill=variable))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = "fill" )+
  scale_x_discrete(labels = c("Samples (n=344)","Blank (n=9)","Blank-seq (n=3)"))+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position="right",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.key.height = unit(1, "mm"))+ 
  guides(fill=guide_legend(ncol=1), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~MetaG[which(MetaG$Matching == "Matched" & MetaG$SC_AMR_alt != "0") ,]$Aetiology_short, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

# Blank_solo<-ggplot(data=TaxaBlank[which(TaxaBlank$Type != "Sample"),],aes(x=Sample_ID, y=value, fill=variable))+
#   scale_fill_manual(values = col_vector_spec) +
#   geom_bar(aes(), stat="identity" )+
#   #scale_y_continuous(labels = scales::percent)+
#   theme(legend.position="right",
#         #axis.text=element_blank(),
#         #axis.title=element_blank(),
#         axis.title=element_text(size=14),
#         #axis.text.x = element_blank(),
#         #axis.text.x = element_text(angle = 90),
#         panel.background = element_rect(fill = NA),
#         axis.line = element_line(size = 0.5, colour = "black"),
#         legend.title = element_blank(),
#         legend.text = element_text(face = "italic"))+ 
#   guides(fill=guide_legend(ncol=1), size = .1)+
#   xlab("")+
#   ylab("Relative abundance (%)")+
#   #facet_wrap(~MetaG[which(MetaG$Matching == "Matched" & MetaG$SC_AMR_alt != "0") ,]$Aetiology_short, scales="free_x")+
#   theme(
#     strip.background = element_rect(
#       color="white", fill="white", size=1, linetype="solid"),
#     strip.text.x = element_text(size = 12)
#   )
# 
# Blank_amr<-ggplot(data=amrBlank,aes(x=Type, y=value, fill=variable))+
#   scale_fill_manual(values = col_vector_spec) +
#   geom_bar(aes(), stat="identity", position = "fill" )+
#   scale_y_continuous(labels = scales::percent)+
#   theme(legend.position="right",
#         #axis.text=element_blank(),
#         #axis.title=element_blank(),
#         axis.title=element_text(size=14),
#         #axis.text.x = element_blank(),
#         #axis.text.x = element_text(angle = 90),
#         panel.background = element_rect(fill = NA),
#         axis.line = element_line(size = 0.5, colour = "black"),
#         legend.title = element_blank(),
#         legend.text = element_text(face = "italic"))+ 
#   guides(fill=guide_legend(ncol=1), size = .1)+
#   xlab("")+
#   ylab("Relative abundance (%)")+
#   #facet_wrap(~MetaG[which(MetaG$Matching == "Matched" & MetaG$SC_AMR_alt != "0") ,]$Aetiology_short, scales="free_x")+
#   theme(
#     strip.background = element_rect(
#       color="white", fill="white", size=1, linetype="solid"),
#     strip.text.x = element_text(size = 12)
#   )

#among 9 blank NC samples only 3 AMR genes were identified pmrA (FQ), tet_B (TC), BRO-2 (BL) in 2 samples. 

pdf(file = "./R_output_files/Figures/Figure_S1.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6)
Blank_taxa 
dev.off()


#Figure 2####
##1. Geography_Resistome####
Geography_Resistome<-ggplot(data=AMRFam,aes(x=Continent, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position="none",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~BphageFamMT$Continent, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##2. Antibiotics_Resistome####
Antibiotics_Resistome<-ggplot(data=AMRFam,aes(x=as.factor(Long.term.antibiotics), y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('No','Yes'))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~BphageFamMT$Continent, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##3. Exacerbation_Resistome###
Exacerbation_Resistome<-ggplot(data=AMRFam[which(is.na(AMRFam$Severity) == FALSE),],aes(x=ExacerbatorState, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('N.Ex (0)','Ex (1-2)', 'F.Ex (3+)'))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~BphageFamMT$Continent, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##4. Severity_Resistome####
Severity_Resistome<-ggplot(data=AMRFam[which(is.na(AMRFam$Severity) == FALSE),],aes(x=Severity, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position="right",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~BphageFamMT$Continent, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##5. Aetiology_Resistome ####
Aetiology_Resistome_Fig3<-ggplot(data=AMRFam[which(is.na(AMRFam$Severity) == FALSE),],aes(x=Aetiology_short, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_grid(~AMRFam$Country, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##6. Microbiology_Resistome ####
Microbiology_Resistome<-ggplot(data=AMRFam,aes(x=Aetiology_short, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('Pseu.','Haem.', 'Strep.', 'other'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_wrap(~AMRFam$PseuDOM, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##7. FEV1 ####
AMRFamMTfev<-subset(AMRFam, FEVfactor != "NA")
FEV1amr<-ggplot(data=AMRFamMTfev,aes(x=FEVfactor, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('>70%','50-70%','30-50%','<30%'))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##8. Virome_AMR####
PhageType_Resistome<-ggplot(data=AMRFam,aes(x=as.factor(SC_AMR_alt), y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('R0','R1', 'R2'))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(ncol=1), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~BphageFamMT$Continent, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

PhageType_Resistome2<-ggplot(data=BphageFam,aes(x=as.factor(SC_AMR_alt), y=RPKM, fill=Virome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#EBA5F3","#5CA5DB","#a3d9d2","#97809e","red","#C6DFA6","#e3ce81","#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#fc5017","#db6960","#B60004","#91CE59","#FFBC06","#3B3B3B","#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('R0','R1', 'R2'))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(ncol=1), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~BphageFamMT$Continent, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##9. Combine resistome ####
#get_legend(Exacerbation_Resistome)->common.leg
Figure_2<-ggarrange(#Geography_Resistome, 
  Exacerbation_Resistome, 
  FEV1amr,
  Severity_Resistome, 
  Aetiology_Resistome_Fig3, 
  #Microbiology_Resistome, 
  #Antibiotics_Resistome, 
  font.label = list(size = 5),
  common.legend = TRUE, legend = "bottom", nrow = 1 ) 
pdf(file = "R_output_files/Figures/Figure_2.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4)
Figure_2
dev.off()

#Figure 3####
AMR_diversity <- Master %>%
  as_tibble() %>%
  select(1:1,64:314) #for genes
  #select(1:1,43:63) #for amr drug class
NAMES_list <- AMR_diversity$SampleID
main_data <- AMR_diversity[AMR_diversity$SampleID %in% NAMES_list, ]
AMR_diversity<-as.matrix(AMR_diversity)
rownames(AMR_diversity) <- AMR_diversity[,1]
AMR_diversity = as.data.frame(subset(AMR_diversity, select = -c(SampleID) ))
AMR_diversity[] <- lapply(AMR_diversity, as.numeric)
AMR_diversity<-AMR_diversity[row.names(AMR_diversity) != "TBS672", , drop = FALSE]

isZero <- base::rowSums(AMR_diversity) == 0
sum(isZero)
AMR_diversity<-AMR_diversity[!isZero,]

MasterVIZ = Master
MasterVIZ$select <- ifelse(MasterVIZ$SC_AMR_alt==0, "null", "Bronchiectasis")
MasterVIZ$select <- ifelse(is.na(MasterVIZ$select), "Non-diseased", MasterVIZ$select)
MasterVIZ$SC_AMR_alt <- ifelse(is.na(MasterVIZ$SC_AMR_alt), "Non-diseased", MasterVIZ$SC_AMR_alt)
AMRDiversityViz<-subset(MasterVIZ, select != "null")
AMRDiversityViz<-AMRDiversityViz[AMRDiversityViz$SampleID != "TBS153", , drop = FALSE] #remove for gene level analysis

AMRDiversityViz_Geo<-subset(AMRDiversityViz, is.na(AMRDiversityViz$SC_AMR_alt) == FALSE & AMRDiversityViz$Matching == "Matched" )

AMR_diversity <- AMR_diversity[ row.names(AMR_diversity) %in% AMRDiversityViz_Geo$SampleID, ]

vegdist(AMR_diversity, "bray")-> Mbiome_PCoA
as.matrix(Mbiome_PCoA)->Mbiome_PCoA
BrayCurtMbiome=cmdscale(Mbiome_PCoA)
#ordiplot (BrayCurtMbiome, display = 'species', type = 'text')
BCords<-scores(BrayCurtMbiome)
BCords<-(as.data.frame(t(BCords)))
BCords<-as.data.frame(t(BCords))

AMRDiversityViz_Geo$Dim1<-BCords$Dim1
AMRDiversityViz_Geo$Dim2<-BCords$Dim2

AMRDiversityViz_Geo$Country <- factor(AMRDiversityViz_Geo$Country, levels = c("SG", "KL", "DD", "MI"))
AMRDiversityViz_Geo$Aetiology_short<- factor(AMRDiversityViz_Geo$Aetiology_short, levels=c("idiopathic", "postInfect", "postTB", "other"))

##AMR PCOA of Resistotypes BY Country####
gg <- data.frame(cluster=factor(AMRDiversityViz_Geo$Country), x=AMRDiversityViz_Geo$Dim1, y=AMRDiversityViz_Geo$Dim2, grp=AMRDiversityViz_Geo$Country)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
PCA_Geo<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2, alpha = 0.5) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_colour_manual(labels = c("Singapore", "Kuala Lumpur", "Dundee", "Milan"), values = c("#CD2C1E","#F7CD46","#2A64B7","#91C55A" ))+
  labs(colour="",  
       x = "PC 1 (23.8%)", y = "PC 2 (4.5%)")+
  theme(legend.position=c(0.8,0.3),
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_y_reverse()+
  xlab("")+
  ylab("")

#is aetiology sigificant?
adonis2(AMR_diversity ~ Aetiology_short, data=AMRDiversityViz_Geo, method="bray", permutations=999)

#Does continent mater more than aetiology?
adonis2(AMR_diversity ~ Continent+Aetiology_short, data=AMRDiversityViz_Geo, method="bray", permutations=999)

#Does country mater more than aetiology?
adonis2(AMR_diversity ~ Country+Aetiology_short, data=AMRDiversityViz_Geo, method="bray", permutations=999)

#Does aetiology mater when assing at gene level?
adonis2(AMR_diversity ~ Aetiology_short, data=AMRDiversityViz_Geo, method="bray", permutations=999)

adonis2(AMR_diversity ~ Continent, data=AMRDiversityViz_Geo, method="bray", permutations=999)


##AMR PCOA of Resistotypes BY Aetiology####
gg <- data.frame(cluster=factor(AMRDiversityViz_Geo$Aetiology_short), x=AMRDiversityViz_Geo$Dim1, y=AMRDiversityViz_Geo$Dim2, grp=AMRDiversityViz_Geo$Aetiology_short)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
PCA_Aet<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2, alpha = 0.5) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_colour_manual(labels = c("IP", "PI", "PTB", "Other"), values = c("#8300c4","#F26B38","#2F9599","grey" ))+
  labs(colour="",  
       x = "PC 1 (23.8%)", y = "PC 2 (4.5%)")+
  theme(legend.position=c(0.7,0.3),
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_y_reverse()+
  xlab("")+
  ylab("")

adonis2(AMR_diversity ~ Aetiology_short, data=AMRDiversityViz_Geo, method="bray", permutations=999)
adonis2(AMR_diversity ~ Country, data=AMRDiversityViz_Geo, method="bray", permutations=999)
adonis2(AMR_diversity ~ Continent, data=AMRDiversityViz_Geo, method="bray", permutations=999)

##Plot the Resistome data####
AMR_geo<-ggplot(data=AMRFam[which(AMRFam$Matching == "Matched") ,],aes(x=Country, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=6), size = .1)+
  xlab("")+
  ylab("")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

AMR_Aet<-ggplot(data=AMRFam[which(AMRFam$Matching == "Matched") ,],aes(x=Aetiology_short, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=6), size = .1)+
  xlab("")+
  ylab("")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##Plot the Taxonomy data####
#wrangle Metagenomic Taxonomy data
MetaG<-read.csv("R_input_files/CAMEB2_bacteria_top50.csv") %>%
  as_tibble() #%>%   
#filter(ExacerbatorState != "NA") %>%
#filter(Matching == "Matched" )

MetaG$SC_AMR_alt<-Master[which(Master$Matching != "NA"),]$SC_AMR_alt
MetaG$Reads<-Master[which(Master$Matching != "NA"),]$ReadsNonHuman
MetaG$Aetiology_short<-Master[which(Master$Matching != "NA"),]$Aetiology_short
MetaG$SampleID <- factor(MetaG$SampleID, levels = MetaG$SampleID[order(MetaG$Reads)])
MetaG$Country<-Master[which(Master$Matching != "NA"),]$Country
MetaG$Country <- factor(MetaG$Country, levels=c("SG", "KL", "DD", "MI"))
MetaG$Aetiology_short <- factor(MetaG$Aetiology_short, levels =c("idiopathic", "postInfect", "postTB", "other"))

MetaG<-melt(MetaG, id.vars = c("SampleID", "Country",  "Continent", "Matching","ExacerbatorState",  "FEV1",   "BSI", "ICS.use", "Oral.ab",   "BMI", "Aetiology", "MMRC.score" ,"DomSpAll", "DomBactOth","DomBactAll", "Reads", "Aetiology_short", "SC_AMR_alt"))

n <- 41
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_spec<-replace(col_vector, 41, "grey") #

Taxa_Aet<-ggplot(data=MetaG[which(MetaG$Matching == "Matched") ,],aes(x=Aetiology_short, y=value, fill=variable))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = "fill" )+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=6), size = .1)+
  xlab("")+
  ylab("")+
  #facet_wrap(~MetaG[which(MetaG$Matching == "Matched" & MetaG$SC_AMR_alt != "0") ,]$Aetiology_short, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

Taxa_Geo<-ggplot(data=MetaG,aes(x=Country, y=value, fill=variable))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = "fill" )+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=6), size = .1)+
  xlab("")+
  ylab("")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

##Plot the phage data####
Phage_Aet<-ggplot(data=BphageFam[which(BphageFam$Matching == "Matched"),],aes(x=Aetiology_short, y=RPKM, fill=Virome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#EBA5F3","#5CA5DB","#a3d9d2","#97809e","red","#C6DFA6","#e3ce81","#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#fc5017","#db6960","#B60004","#91CE59","#FFBC06","#3B3B3B","#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

Phage_Geo<-ggplot(data=BphageFam[which(BphageFam$Matching == "Matched"),],aes(x=Country, y=RPKM, fill=Virome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#EBA5F3","#5CA5DB","#a3d9d2","#97809e","red","#C6DFA6","#e3ce81","#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#fc5017","#db6960","#B60004","#91CE59","#FFBC06","#3B3B3B","#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=2), size = .1)+
  xlab("")+
  ylab("")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

Figure_3<-ggarrange(PCA_Geo,Taxa_Geo,Phage_Geo,
                    PCA_Aet,Taxa_Aet,Phage_Aet,
  font.label = list(size = 5),
  common.legend = FALSE, ncol = 3,nrow = 2)

Figure_3_Legend_CF<-Taxa_Aet+
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        legend.key.height = unit(1, "mm"))+
  guides(fill=guide_legend(ncol=1), size = .1)

Figure_3_Legend_DG<-Phage_Aet+theme(legend.position = "right")+
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        legend.key.height = unit(1, "mm"))+
  guides(fill=guide_legend(ncol=1), size = .1)

pdf(file = "R_output_files/Figures/Figure_3.pdf",   # The directory you want to save the file in
   width = 10, # The width of the plot in inches
  height = 6)
Figure_3 
dev.off()

pdf(file = "./R_output_files/Figures/Figure_3_Legend_CF.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
Figure_3_Legend_CF 
dev.off()

pdf(file = "./R_output_files/Figures/Figure_3_Legend_DG.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
Figure_3_Legend_DG
dev.off()





#Figure 4#####
#AMR PCOA of Resistotypes BY SC_RESISTOTYPE   
##Panel A####
##Panel B####
#AMRDiversityViz_Geo<-as.data.frame(merge(AMRDiversityViz_Geo, lab, by.x = 1, by.y = 0, all.x = TRUE)) #add in labels
gg <- data.frame(cluster=factor(AMRDiversityViz_Geo$SC_AMR_alt), x=AMRDiversityViz_Geo$Dim1, y=AMRDiversityViz_Geo$Dim2, grp=AMRDiversityViz_Geo$SC_AMR_alt)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
PCA_RT<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2, alpha = 0.5) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_colour_manual(values = c("#1800F5","#932DE7"), labels = c("RT1", "RT2"))+
  labs(colour="",  
       x = "PC 1 (23.8%)", y = "PC 2 (4.5%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_y_reverse()+
  #scale_x_reverse()+
  guides(colour = guide_legend(reverse = FALSE))

##Panel C####

##Panel DEF####
RT1_RT2_plot<-AMRDiversityViz %>%
  subset(select == "Bronchiectasis") %>%
  subset(Matching == "Matched")

Ex<-ggplot(data=RT1_RT2_plot, aes(as.factor(SC_AMR_alt), Exacerbations, group = as.factor(SC_AMR_alt), fill = as.factor(SC_AMR_alt)))+ 
  scale_fill_manual(values=c("#1800F5","#932DE7"))+ 
  scale_y_continuous(breaks=seq(0,12,2))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0.1)+
  facet_grid(. ~ "Exacerbations")+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect( fill="white",size = 1))+
  labs(fill='Resistotype')+
  scale_x_discrete(labels = c('RT1', 'RT2'))

FEV1<-ggplot(data=RT1_RT2_plot, aes(as.factor(SC_AMR_alt), FEV1, group = as.factor(SC_AMR_alt), fill = as.factor(SC_AMR_alt)))+ 
  scale_fill_manual(values=c("#1800F5","#932DE7"))+ 
  expand_limits(y=c(0,175))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0.1)+
  facet_wrap(~ "FEV1 (% predicted)")+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect( fill="white",size = 1))+
  labs(fill='Resistotype')+
  scale_x_discrete(labels = c('RT1', 'RT2'))

BSI<-ggplot(data=RT1_RT2_plot, aes(as.factor(SC_AMR_alt), BSI, group = as.factor(SC_AMR_alt), fill = as.factor(SC_AMR_alt)))+ 
  scale_fill_manual(values=c("#1800F5","#932DE7"))+ 
  expand_limits(y=c(0,22))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0.1)+
  facet_grid(. ~ "Severity (BSI)")+
  theme(legend.position="none",
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.title = element_blank(),
    axis.line = element_line(size = 0.5, colour = "black"),
    panel.background = element_rect(fill = NA),
    strip.background = element_rect(fill="white",size = 1))+
  labs(fill='Resistotype')+
  scale_x_discrete(labels = c('RT1', 'RT2'))

Figure_4DEF<-ggarrange(Ex,FEV1, BSI, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

pdf(file = "./R_output_files/Figures/Figure_4DEF.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
Figure_4DEF
dev.off()

pdf(file = "./R_output_files/Figures/Figure_4B.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
PCA_RT
dev.off()



#some statistical analysis [required?]
# wilcox.test(RT1_RT2_plot[which(RT1_RT2_plot$Continent == "Europe"),]$FEV1~RT1_RT2_plot[which(RT1_RT2_plot$Continent == "Europe"),]$SC_AMR_alt)
# wilcox.test(RT1_RT2_plot$Exacerbations~RT1_RT2_plot$SC_AMR_alt)
# wilcox.test(RT1_RT2_plot$FEV1~RT1_RT2_plot$SC_AMR_alt)
# wilcox.test(RT1_RT2_plot$BSI~RT1_RT2_plot$SC_AMR_alt)
# wilcox.test(RT1_RT2_plot$Age~RT1_RT2_plot$SC_AMR_alt)
# chisq.test(RT1_RT2_plot$Continent,RT1_RT2_plot$SC_AMR_alt)
# chisq.test(RT1_RT2_plot$Country,RT1_RT2_plot$SC_AMR_alt)
# chisq.test(RT1_RT2_plot$Country,RT1_RT2_plot$Long.term.antibiotics)
# chisq.test(RT1_RT2_plot$Country,RT1_RT2_plot$PseuDOM)
# chisq.test(RT1_RT2_plot$Continent,RT1_RT2_plot$Aetiology_short)
# chisq.test(RT1_RT2_plot$Country,RT1_RT2_plot$Aetiology_short)
# table(RT1_RT2_plot$DomBactAll,RT1_RT2_plot$Continent)
# chisq.test(MasterMT$Continent,MasterMT$Aetiology_short)
# summary(MasterMT$Exacerbations,MasterMT$Continent)
# wilcox.test(MasterMT$Sex..Male.0..Female.1.~MasterMT$Continent)
# summary(MasterMT[which(MasterMT$Continent == "Europe"),]$Severity,
#         MasterMT[which(MasterMT$Continent == "Europe"),]$Continent)
# table(MasterMT$Severity,
#         MasterMT$Aetiology_short)
# table(MasterMT[which(MasterMT$Continent == "Asia"),]$Continent,
#       MasterMT[which(MasterMT$Continent == "Asia"),]$Aetiology_short)
# summary(data=MasterMT[which, MasterMT$Continent == "Europe"],MasterMT$Exacerbations)
# table(RT1_RT2_plot$Country,RT1_RT2_plot$PseuDOM)
# table(RT1_RT2_plot$Continent,RT1_RT2_plot$Aetiology_short)
# table(RT1_RT2_plot$Country,RT1_RT2_plot$DomBactOth_short)
# table(RT1_RT2_plot$SC_AMR_alt,RT1_RT2_plot$Number.of.lobes.affected)

#Seems belowis and bar plot for AMR * country [required ?]
# RT_AMR<-ggplot(data=AMRFam[which(AMRFam$CTRL=="PATIENT" & AMRFam$SC_AMR_alt != "0"),],aes(x=Country, y=RPKM, fill=Resistome))+
#   geom_bar(aes(), stat="identity", position = "fill") +
#   scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
#   scale_y_continuous(labels = scales::percent)+
#   #scale_x_discrete(labels = c('RT1','RT2'))+
#   theme(legend.position="bottom",
#         #axis.text=element_blank(),
#         #axis.title=element_blank(),
#         axis.title=element_text(size=14),
#         #axis.text.x = element_blank(),
#         #axis.text.x = element_text(angle = 90),
#         panel.background = element_rect(fill = NA),
#         axis.line = element_line(size = 0.5, colour = "black"),
#         legend.title = element_blank(),
#         legend.text = element_text(face = "italic"))+ 
#   guides(fill=guide_legend(nrow=3), size = .1)+
#   xlab("")+
#   ylab("Relative abundance (%)")+
#   facet_wrap(~AMRFam[which(AMRFam$CTRL=="PATIENT" & AMRFam$SC_AMR_alt != "0"),]$CTRL, scales="free_x")+
#   theme(
#     strip.background = element_rect(
#       color="white", fill="white", size=1, linetype="solid"),
#     strip.text.x = element_text(size = 12)
#   )
# MetaG$Aetiology_short<- factor(MetaG$Aetiology_short, levels=c("idiopathic", "postInfect", "postTB", "other"))

#Figure 5 ####
Taxa_Geo_Fig5a<-ggplot(data=MetaG[which(MetaG$SC_AMR_alt != "0"),],aes(x=as.factor(SC_AMR_alt), y=value, fill=variable))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = "fill" )+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('RT1','RT2'))+
  ggtitle("CAMEB2 (n=251)") +
  theme(
    legend.position="none",
    #axis.text=element_blank(),
    #axis.title=element_blank(),
    axis.title=element_text(size=14),
    #axis.text.x = element_blank(),
    #axis.text.x = element_text(angle = 90),
    panel.background = element_rect(fill = NA),
    axis.line = element_line(size = 0.5, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(face = "italic", size = 7),
    legend.key.height = unit(1, "mm"))+ 
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~MetaG[which(MetaG$SC_AMR_alt != "0"),]$SC_AMR_alt, scales="free_x",nrow = 1)+
  theme(plot.title = element_text(hjust = 0.5, size = 12), 
        plot.title.position = "plot", 
        #plot.title.margin = margin(b = 10),
        strip.background = element_rect( 
          color="white", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12)
  )+
  guides(fill=guide_legend(ncol=1), size = .1)

variable_names_a <- list(
  "Asia" = "SG-KL (n=130)",
  "Europe" = "DD-MI (n=121)"
)

variable_labeller <- function(variable,value){
  return(variable_names_a[value])
}

Taxa_Geo_Fig5b<-ggplot(data=MetaG[which(MetaG$SC_AMR_alt != "0"),],aes(x=as.factor(SC_AMR_alt), y=value, fill=variable))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = "fill" )+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('RT1','RT2'))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 7),
        legend.key.height = unit(1, "mm"))+ 
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_wrap(~MetaG[which(MetaG$SC_AMR_alt != "0"),]$Continent, scales="free_x",nrow = 1,labeller=variable_labeller)+
  theme(
    strip.background = element_rect( 
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )+
  guides(fill=guide_legend(ncol=1), size = .1)

variable_names_b<- list(
  "SG" = "SG (n=95)",
  "KL" = "KL (n=35)",
  "DD" = "DD (n=96)",
  "MI" = "MI (n=25)"
)

variable_labeller <- function(variable,value){
  return(variable_names_b[value])
}

Taxa_Geo_Fig5c<-ggplot(data=MetaG[which(MetaG$SC_AMR_alt != "0"),],aes(x=as.factor(SC_AMR_alt), y=value, fill=variable))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = "fill" )+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c('RT1','RT2'))+
  theme(legend.position="none",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 7),
        legend.key.height = unit(1, "mm"))+ 
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_wrap(~MetaG[which(MetaG$SC_AMR_alt != "0"),]$Country, scales="free_x",nrow = 1,labeller=variable_labeller)+
  theme(
    strip.background = element_rect( 
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )+
  guides(fill=guide_legend(ncol=1), size = .1)


# Combine the first two plots horizontally
plot12 <- ggarrange(Taxa_Geo_Fig5a, Taxa_Geo_Fig5b, ncol = 2, widths = c(1, 1))
Figure_5ABC<-ggarrange(plot12, Taxa_Geo_Fig5c, ncol = 1, heights = c(1, 1), common.legend = TRUE, legend = 'right')

## BIPARTITE GRAPHS ####
### Species-AMR ####
#### Project name and working directory ####
project_name <- "CAMEB2"
wkdir <- file.path("./R_input_files/Ivan_AMR_analysis_2020", project_name)

indir <- file.path(wkdir, "analysis", "AMR-CONTIGS", "unweighted-bacteria-phages-plasmids")
outdir <- indir

# AMR-contigs
df <- read.csv(file.path(indir, "AMR-contigs-Species-abundance.csv"), row.names= 1, check.names= FALSE)
df <- df[df$Type=="Chromosome", ]
df <- df[df$Species!="Unknown sp.", ]

SpeciesAMR <- data.frame(Species= gsub(" - .*", "", unique(df$SpeciesAMR)), AMR= gsub(".* - ", "", unique(df$SpeciesAMR)))
AMRcountsPerSpecies <- data.frame(sort(table(SpeciesAMR$Species))); colnames(AMRcountsPerSpecies) <- c("Species", "Count")
SpeciesCountsPerAMR <- data.frame(sort(table(SpeciesAMR$AMR))); colnames(SpeciesCountsPerAMR) <- c("AMR", "Count")

kaiju<-readRDS("./R_input_files/Ivan_AMR_analysis_2020/CAMEB2/analysis/KAIJU/ps_Species.RData")

ps.prop <- transform_sample_counts(kaiju, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
otutable<-as.data.frame(otu_table(ps.prop))
z=colSums(otutable>0.1) #filter to reduced number of taxa 
sel_col=row.names(as.data.frame(z[z>=(0.1*(nrow(otutable)))])) #In 5% patients prevalent
otutable<-otutable[sel_col]
remove(sel_col,z)
otutable["Sample_ID"]<-row.names(otutable)
otutable<-gather(otutable, Bacteriome, value, -Sample_ID)
taxa.means<-aggregate(value~Bacteriome, FUN=mean, data =otutable)
taxa.means<-taxa.means[
  with(taxa.means, order(-value)),
]
taxa.order<-as.vector(taxa.means$Bacteriome)
taxa.order.top50<-taxa.order[1:42]

#write.csv(otutable, file.path(wkdir, 'SP_otutable.csv'))
#write.csv(SpeciesAMR, file.path(outdir, "SpeciesAMR.csv"), row.names= FALSE)
#write.csv(AMRcountsPerSpecies, file.path(outdir, "SpeciesAMR-AMRcountsPerSpecies.csv"))
#write.csv(SpeciesCountsPerAMR, file.path(outdir, "SpeciesAMR-SpeciesCountsPerAMR.csv"))

# Barplot SpeciesAMR-AMRcountsPerSpecies
# ----------------------------------------------------------------------------------------------

#need to pull taxa names from master datafile for CAMEB2. Obviously there must be a better way.
Master <-read.csv("./R_input_files/MetaDataWsnfResults_3.4_TJ.csv") %>% # v 3.3 includes raw read data
  as_tibble()

## Read files
colTaxa <- readRDS(file.path(wkdir, "analysis", "KAIJU", "colTaxa_Species.RData"))
names(colTaxa) <- gsub("_", " ", names(colTaxa))

##alt colours [original CAMEB2 colours]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_spec<-replace(col_vector, 44, "grey")
cols_AMR<-col_vector[1:20]
names(cols_AMR) <- gsub("\\.", " ", colnames(Master[646:665]))

data <- read.csv(file.path(indir, "SpeciesAMR-AMRcountsPerSpecies.csv"))
data <- data[data$Species %in% as.list(taxa.order.top50), ]
data <- merge(data, taxa.means, by.x = "Species", 
              by.y = "Bacteriome", all.x = TRUE, all.y = FALSE)
data$Species <- factor(data$Species, levels = data$Species[order(data$Count)])
colSpecies <- append((rep("gray40", nrow(data)-20)),rev(col_vector[1:20]))
names(colSpecies) <- as.vector(data$Species[order(data$value)])
colSpecies<-colSpecies[order(match(names(colSpecies),data$Species[order(data$Count)]))]



## Plotting
p <- ggplot(data, aes(x=Species, y=Count, fill=Species)) + 
  geom_bar(stat = "identity") + #ggtitle("Number of AMR genes associated with bacterial species") +
  coord_flip() + xlab("") + ylab("N AMR genes") +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color="black", size= 6, angle= 0, hjust= 1),
        axis.text.y = element_text(face="bold.italic", size= 8, colour=colSpecies),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color="gray4", fill=NA),
        #panel.grid.major = element_line(color="gray80"),
        panel.grid.minor = element_line(color="gray80")
  ) + scale_x_discrete(position = "bottom") +
  scale_fill_manual(name= "Species", values=colSpecies)
p

#ggsave(file=file.path(outdir, "SpeciesAMR-AMRcountsPerSpecies.png"), width = 120, height = 300, units = "mm", dpi = 500)

#-------------------

shortbread<-readRDS("./R_input_files/Ivan_AMR_analysis_2020/CAMEB2/analysis/AMR/ps_Gene.RData")
ps.prop <- transform_sample_counts(shortbread, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
otutable<-as.data.frame(otu_table(ps.prop))
z=colSums(otutable>0.1) #filter to reduced number of taxa 
sel_col=row.names(as.data.frame(z[z>=(0.1*(nrow(otutable)))])) #In 5% patients prevalent
otutable<-otutable[sel_col]
remove(sel_col,z)
otutable["Sample_ID"]<-row.names(otutable)
otutable<-gather(otutable, Resistome, value, -Sample_ID)
taxa.means<-aggregate(value~Resistome, FUN=mean, data =otutable)
taxa.means<-taxa.means[
  with(taxa.means, order(-value)),
]
ps <- readRDS(file.path(wkdir, "analysis", "AMR", "ps_Gene.RData"))
taxa.means$Resistome <- c(paste0(tax_table(ps)[setdiff(taxa.means$Resistome, "Others"),"Gene"], " (", tax_table(ps)[setdiff(taxa.means$Resistome, "Others"),"Drug"], ")"))
taxa.order<-as.vector(taxa.means$Resistome)
taxa.order.top50<-taxa.order[1:50] #40 genes in final graph

# Barplot SpeciesAMR-SpeciesCountsPerAMR
# ----------------------------------------------------------------------------------------------

colTaxa <- readRDS(file.path(wkdir, "analysis", "AMR", "colTaxa_Gene.RData"))
#names(colTaxa) <- c(tax_table(ps)[setdiff(names(colTaxa), "Others"),"Gene"])
names(colTaxa) <- c(paste0(tax_table(ps)[setdiff(names(colTaxa), "Others"),"Gene"], " (", tax_table(ps)[setdiff(names(colTaxa), "Others"),"Drug"], ")"))

data <- read.csv(file.path(indir, "SpeciesAMR-SpeciesCountsPerAMR.csv"))
#data$AMR <- factor(gsub(" \\(.*", "", as.vector(data$AMR)), levels= gsub(" \\(.*", "", as.vector(data$AMR)))
data$AMR <- factor(as.vector(data$AMR), levels= as.vector(data$AMR))
data <- data[data$AMR %in% as.list(taxa.order.top50), ]
data <- merge(data, taxa.means, by.x = "AMR", 
              by.y = "Resistome", all.x = TRUE, all.y = FALSE)
#data$AMR <- factor(data$AMR, levels = data$AMR[order(data$Count)])
colAMR <- rep("gray40", nrow(data))
names(colAMR) <- as.vector(data$AMR)
colAMR[names(colTaxa)] <- colTaxa
colAMR<-colAMR[order(match(names(colAMR),data$AMR[order(data$Count)]))]

#data<-data[!grepl("Streptomyces rishiriensis parY", data$AMR),]

p2 <- ggplot(data, aes(x=AMR, y=Count, fill=AMR)) + 
  geom_bar(stat = "identity") + #ggtitle("Number of bacterial species associated with AMR genes") +
  coord_flip() + xlab("") + ylab("N bacteria") +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color="black", size= 6, angle= 0, hjust= 1),
        axis.text.y = element_text(face="bold.italic", size= 8, colour=colAMR),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color="gray4", fill=NA),
        #panel.grid.major = element_line(color="gray80"),
        panel.grid.minor = element_line(color="gray80")
  ) +
  scale_fill_manual(name= "AMR", values=colAMR)

p2

#ggsave(file=file.path(outdir, "SpeciesAMR-SpeciesCountsPerAMR.png"), width = 183, height = 300, units = "mm", dpi = 500)

pTEMP<-ggarrange(p,p2, common.legend = FALSE)
Figure_5<-ggarrange(Figure_5ABC,pTEMP, ncol = 1, heights = c(1, 0.8),common.legend = FALSE)
pdf(file = "./R_output_files/Figures/Figure_5.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 14)
Figure_5
dev.off()

#Figure S3#####
##ICS_stack####
ICSanalysis<-ggplot(data=AMRFam[which(is.na(AMRFam$Matching) == FALSE) ,],aes(x=as.factor(ICS.use), y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c("No ICS","ICS"))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=3), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )+
  ggtitle("Inhaled corticosteroid use")

##ABX_stack####
ABXanalysis<-ggplot(data=AMRFam[which(is.na(AMRFam$Matching) == FALSE) ,],aes(x=as.factor(Long.term.antibiotics
), y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c("No macrolide","macrolide"))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=3), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_wrap(~AMRFam$CTRL, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )+
  ggtitle("Long-term macrolide use")

## PCA plot ####
cohort<-subset(Master, is.na(Master$SC_AMR_alt) == FALSE & SC_AMR_alt != 0 & Matching == "Matched" )

AMR_diversity <- cohort %>%
  as_tibble() %>%
  #select(1:1,395:645) #for genes
  select(1:1,374:394) #for amr drug class
NAMES_list <- AMR_diversity$SampleID
main_data <- AMR_diversity[AMR_diversity$SampleID %in% NAMES_list, ]
AMR_diversity<-as.matrix(AMR_diversity)
rownames(AMR_diversity) <- AMR_diversity[,1]
AMR_diversity = as.data.frame(subset(AMR_diversity, select = -c(SampleID) ))
AMR_diversity[] <- lapply(AMR_diversity, as.numeric)
AMR_diversity<-AMR_diversity[row.names(AMR_diversity) != "TBS672", , drop = FALSE]

isZero <- base::rowSums(AMR_diversity) == 0
sum(isZero)
AMR_diversity<-AMR_diversity[!isZero,]

vegdist(AMR_diversity, "bray")-> Mbiome_PCoA
as.matrix(Mbiome_PCoA)->Mbiome_PCoA
BrayCurtMbiome=cmdscale(Mbiome_PCoA)
ordiplot (BrayCurtMbiome, display = 'species', type = 'text')
BCords<-scores(BrayCurtMbiome)
BCords<-(as.data.frame(t(BCords)))
BCords<-as.data.frame(t(BCords))

MasterVIZ = Master
MasterVIZ$select <- ifelse(MasterVIZ$SC_AMR_alt==0, "null", "Bronchiectasis")
MasterVIZ$select <- ifelse(is.na(MasterVIZ$select), "Non-diseased", MasterVIZ$select)
MasterVIZ$SC_AMR_alt <- ifelse(is.na(MasterVIZ$SC_AMR_alt), "Non-diseased", MasterVIZ$SC_AMR_alt)
AMRDiversityViz<-subset(MasterVIZ, select != "null")
AMRDiversityViz$select

AMRDiversityViz<-AMRDiversityViz[AMRDiversityViz$SampleID != "TBS153", , drop = FALSE] #remove for gene level analysis

AMRDiversityViz <- AMRDiversityViz[AMRDiversityViz$SampleID %in% cohort$SampleID, ]
AMRDiversityViz <- AMRDiversityViz[AMRDiversityViz$SampleID %in% row.names(AMR_diversity), ]

AMRDiversityViz$Dim1<-BCords$Dim1
AMRDiversityViz$Dim2<-BCords$Dim2

#checking PC %s
rda(X = AMR_diversity, scale = TRUE)
checkEig<-capscale(AMR_diversity ~1)
Eig <-eigenvals(checkEig)
Eig / sum(Eig)

####ICS_PCA####   
gg <- data.frame(cluster=factor(AMRDiversityViz$ICS.use), x=AMRDiversityViz$Dim1, y=AMRDiversityViz$Dim2, grp=AMRDiversityViz$ICS.use)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
ICS.pca<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  #scale_shape_discrete(labels = c("Healthy", "Bronchiectasis"))+
  scale_colour_manual(values = c("#F8766D", "#619CFF"), labels = c("ICS","No ICS"))+
  labs(colour="",  
       x = "PC 1 (77.6%)", y = "PC 2 (6.9%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
        legend.key.size = unit(3, 'mm')
  )+
  scale_x_reverse()+
  #scale_y_reverse()+ #add for gene level analysis
  guides(colour = guide_legend(reverse = T))

adonis2(AMR_diversity ~ ICS.use, data=AMRDiversityViz, method="bray", permutations=9999)

####ABX_PCA####   
gg <- data.frame(cluster=factor(AMRDiversityViz$Long.term.antibiotics), x=AMRDiversityViz$Dim1, y=AMRDiversityViz$Dim2, grp=AMRDiversityViz$Long.term.antibiotics)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
ABX.pca<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  #scale_shape_discrete(labels = c("Healthy", "Bronchiectasis"))+
  scale_colour_manual(values = c("#F8766D", "#619CFF"), labels = c("macrolide","No macrolide"))+
  labs(colour="",  
       x = "PC 1 (77.6%)", y = "PC 2 (6.9%)")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
        legend.key.size = unit(3, 'mm')
  )+
  scale_x_reverse()+
  #scale_y_reverse()+ #add for gene level analysis
  guides(colour = guide_legend(reverse = T))

adonis2(AMR_diversity ~ Long.term.antibiotics, data=AMRDiversityViz, method="bray", permutations=9999)



Figure_S3_temp_1<-ggarrange(ICSanalysis, ABXanalysis,common.legend = TRUE, legend ="bottom")
Figure_S3_temp_2<-ggarrange(ICS.pca, ABX.pca,
common.legend = FALSE)
Figure_S3<-ggarrange(Figure_S3_temp_1, Figure_S3_temp_2,common.legend = FALSE, ncol = 1)

pdf(file = "./R_output_files/Figures/Figure_S3.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12)
Figure_S3
dev.off()



#code for a PCA analysis on the above swap country to Long.term.ab.use
gg <- data.frame(cluster=factor(AMRDiversityViz_Geo$ICS.use), x=AMRDiversityViz_Geo$Dim1, y=AMRDiversityViz_Geo$Dim2, grp=AMRDiversityViz_Geo$Country)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
PCA_Geo<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2, alpha = 0.5) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  #scale_colour_manual(labels = c("Singapore", "Kuala Lumpur", "Dundee", "Milan"), values = c("#CD2C1E","#F7CD46","#2A64B7","#91C55A" ))+
  labs(colour="",
       x = "PC 1 (23.8%)", y = "PC 2 (4.5%)")+
  theme(legend.position=c(0.8,0.3),
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_x_reverse()+
  scale_y_reverse()+
  xlab("")+
  ylab("")

adonis2(AMR_diversity ~ ICS.use, data=AMRDiversityViz_Geo, method="bray", permutations=999)
#some betadispers work : Long, William. (2018). Re: How should I correctly manage PERMANOVA for factors with interactions?. Retrieved from: https://www.researchgate.net/post/How_should_I_correctly_manage_PERMANOVA_for_factors_with_interactions/5bc11473f0fb62111712404e/citation/download. 
#https://www.rdocumentation.org/packages/vegan/versions/2.6-4/topics/permutest.betadisper

mod<-betadisper(vegdist(AMR_diversity, "bray"), as.factor(AMRDiversityViz_Geo$Aetiology_short))
anova(mod)
pmod<-permutest(mod,permutations = 99, pairvise = T)
(mod.HSD <-TukeyHSD(mod))
plot(mod.HSD)
pstat<-permustats(pmod)
densityplot(pstat, scales = list(x=list(relation = "free")))
qqmath(pstat, scales  = list(relation = "free"))

# adonis2(AMRLT_diversity ~ Severity, data=LTDiversityViz, method="bray", permutations=999)

pdf(file = "./R_output_files/Figures/Figure_S4.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
ABXanalysis
dev.off()

#Figure 6####
MasterANDGK <-read.csv("../LONGITUDINAL/OneDrive_1_24-06-2021/MetaDataWsnfResults_GK_final.csv") # v 2.0 includes an aggregated 'others' column. v 3.0 amalgamates all data views

cohort<-subset(MasterANDGK, is.na(SC_AMR_alt) == FALSE & SC_AMR_alt !=0)
cohort<-tail(cohort, n =8)
cohort <- cohort %>% #clinical variables + amr families
  as_tibble() %>%
  select(-48:-373,-395:-686)
cohort$PseuDOM<-ifelse(cohort$DomBactOth_short == "Pseudomonas", "Pseudomonas", "other")
cohort$HaemDOM<-ifelse(cohort$DomBactOth_short == "Haemophilus", "Haemophilus", "other")
cohort$FEVfactor<-cut(cohort$FEV1,  breaks=c(0, 30, 50, 70, Inf))
cohort$AetiologyN3<-ifelse(cohort$Aetiology_short == "postTB", "postInfect", cohort$Aetiology_short)
#set levels
cohort$ExacerbatorState <- factor(cohort$ExacerbatorState, levels=c("NonEx", "Exacerbator", "FreqEx"))
cohort$Country <- factor(cohort$Country, levels=c("SG", "KL", "DD", "MI"))
cohort$Aetiology_short <- factor(cohort$Aetiology_short, levels=c("idiopathic", "postInfect", "postTB", "other"))
cohort$DomBactOth_short <- factor(cohort$DomBactOth_short, levels=c("Pseudomonas", "Haemophilus", "Streptococcus", "other"))
cohort$SampleID <- factor(cohort$SampleID, levels = cohort$SampleID[order(cohort$SC_V)])
cohort$FEVfactor<-fct_rev(cohort$FEVfactor)
cohort <- cohort %>% 
  gather(Resistome, RPKM, Acridine.dye,Aminocoumarin.antibiotic,Aminoglycoside,Antibacterial.free.fatty.acids,Beta.lactam,Bicyclomycin,Diaminopyrimidine,Fluoroquinolone,Fosfomycin,Fusidic.acid,MLS,Multidrug,Mupirocin,Nitroimidazole.antibiotic,Nucleoside.antibiotic,Peptide.antibiotic,Phenicol,Rifampicin,Sulfonamide.antibiotic,Tetracycline,Triclosan, -SampleID, -Country,  -Continent, -Matching, -Paired, -Trio, -Age,-Sex..Male.0..Female.1.,-Exacerbations,-ExacerbatorState,  -FEV1,   -BSI, -ICS.use, -Oral.ab,   -BMI, -Aetiology, -Aetiology_short,-MMRC.score, -SC_V, -SC_AMR, -SC_snf_tuned, -SC_VAT_webtool, -SC_AMR_alt, -DomBactOth_short, -PseuDOM, -FEVfactor, -AetiologyN3)
cohort$PseuDOM <- factor(cohort$PseuDOM, levels = c("Pseudomonas", "other"))
cohort$HaemDOM <- factor(cohort$HaemDOM, levels = c("Haemophilus", "other"))
cohort$AetiologyN3 <- factor(cohort$AetiologyN3, levels = c("idiopathic", "postInfect", "other"))
cohortMT<-subset(cohort, Matching == "Matched") #Matched data
cohortPR<-subset(cohort, Paired == "Pair")
cohortTR<-subset(cohort, Trio != "NonTrio")
cohortUM<-subset(cohort, Matching == "UnMatched")
cohort$CTRL<-ifelse(is.na(cohort$Age), "CTRL", "PATIENT")
cohort$ReadsNonHuman <- factor(cohort$ReadsNonHuman, levels=c("Pre","Post"))


AMR_GK<-ggplot(data=cohort,aes(x=ReadsTrimmed, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('Non-diseased','Bronchiectasis'))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 9),
        legend.key.height = unit(1.5, "mm"))+
  guides(fill=guide_legend(ncol=1), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_wrap(~cohort$ReadsNonHuman, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    #strip.text.x = element_blank()
  )

MasterGK <-read.csv("../LONGITUDINAL/OneDrive_1_24-06-2021/GK_master_Taxa25.csv") # drop taxa numbers to 25 + other

AMRGK <- MasterGK %>%
  as_tibble() %>%
  select(1:50)
AMR_cols<-colnames(AMRGK[14:39])
AMR_cols
AMRGK <- AMRGK %>% 
  gather(AMR, RPKM, AMR_cols, -SampleSeqNo, -SputumSampleNo,  -TypeSamples, -TypeSamplesA,-TypeSamplesB,-Exacerbations,	-FEV1, -BSI,	-Severity, -WksToNxtEx,	-TmToNxtEx,	-Antibiotic,	-Antibiotic_class)

AMRGK$AMR <- factor(AMRGK$AMR, levels = AMR_cols )
AMRGK$TypeSamples <- factor(AMRGK$TypeSamples , levels = c("Pre", "Post"))
AMRGK<-AMRGK[is.na(AMRGK$TypeSamples) != TRUE, , drop = FALSE]

n <- 56
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_spec<-replace(col_vector, 26, "grey") #

TAXA_GK<-ggplot(data=AMRGK,aes(x=SputumSampleNo
, y=RPKM, fill=AMR))+
  scale_fill_manual(values = col_vector_spec) +
  geom_bar(aes(), stat="identity", position = 'fill')+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('IP','PI', 'PTB', "other"))+
  theme(legend.position="right",
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 9),
        legend.key.height = unit(1.5, "mm"))+
  guides(fill=guide_legend(ncol=1), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_wrap(~AMRGK$TypeSamples, scales="free_x")+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    #strip.text.x = element_blank()
  )

Master <-read.csv("../LONGITUDINAL/OneDrive_1_24-06-2021/MetaDataWsnfResults_GK_final.csv") %>%
  #Master <-read.csv("R_input_files/MetaDataWsnfResults_3.5temp_TJ.csv") %>% # v 3.3 includes raw read data, 3.4 includes info on paring and "trios"
  as_tibble()

AMR_diversity <- Master %>%
  subset(is.na(SC_AMR_alt)==FALSE) %>%
  as_tibble() %>%
  select(1:1,395:645) #for genes
#select(1:1,374:394) #for amr drug class
NAMES_list <- AMR_diversity$SampleID
main_data <- AMR_diversity[AMR_diversity$SampleID %in% NAMES_list, ]
AMR_diversity<-as.matrix(AMR_diversity)
rownames(AMR_diversity) <- AMR_diversity[,1]
AMR_diversity = as.data.frame(subset(AMR_diversity, select = -c(SampleID) ))
AMR_diversity[] <- lapply(AMR_diversity, as.numeric)
AMR_diversity<-AMR_diversity[row.names(AMR_diversity) != "TBS672", , drop = FALSE]

diversity(AMR_diversity, "shannon")

isZero <- base::rowSums(AMR_diversity) == 0
sum(isZero)
AMR_diversity<-AMR_diversity[!isZero,]

287-37

MasterVIZ = Master
MasterVIZ <- MasterVIZ[ MasterVIZ$SampleID %in% row.names(AMR_diversity), ]

vegdist(AMR_diversity, "bray")-> Mbiome_PCoA
as.matrix(Mbiome_PCoA)->Mbiome_PCoA
BrayCurtMbiome=cmdscale(Mbiome_PCoA)
ordiplot (BrayCurtMbiome, display = 'species', type = 'text')
BCords<-scores(BrayCurtMbiome)
BCords<-(as.data.frame(t(BCords)))
BCords<-as.data.frame(t(BCords))

MasterVIZ$Dim1<-BCords$Dim1
MasterVIZ$Dim2<-BCords$Dim2

MasterVIZ$Country <- factor(MasterVIZ$Country, levels = c("SG", "KL", "DD", "MI"))
MasterVIZ$Aetiology_short<- factor(MasterVIZ$Aetiology_short, levels=c("idiopathic", "postInfect", "postTB", "other"))

MasterVIZ$PrePost<-ifelse(MasterVIZ$ReadsNonHuman == "Pre","Pre",NA)
MasterVIZ$PrePost<-ifelse(MasterVIZ$ReadsNonHuman == "Post","Post", MasterVIZ$PrePost)
MasterVIZ$PrePost<-ifelse(is.na(MasterVIZ$PrePost) == TRUE, MasterVIZ$SC_AMR_alt, MasterVIZ$PrePost)

gg <- data.frame(cluster=factor(MasterVIZ$PrePost), x=MasterVIZ$Dim1, y=MasterVIZ$Dim2, grp=MasterVIZ$SC_AMR_alt)
# calculate group centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# generate star plot...
PA_erd_PCA<-ggplot(gg) +
  #scale_col_manual(values=c(16, 16, 16,16))+
  scale_linetype_identity() +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster),alpha = 0.3)+
  geom_point(aes(x=x,y=y, colour = cluster), size = 2, alpha = 0.5) + #can add ",shape = shape" in aes to introduce shape to points.
  #geom_point(aes(x=x,y=y, colour = cluster, shape = shape), size = 2) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5) +
  geom_point(data=centroids, aes(x=x, y=y, color=cluster), size=5, shape = 13, colour = "black") +
  scale_colour_manual(values = c("#1800F5","#932DE7", "red", "green"), labels = c("RT1", "RT2", "Post-eradication", "Pre-eradication"))+
  labs(colour="",  
       x = "PC 1 (23.8%)", y = "PC 2 (4.5%)")+
  theme(legend.position=c(0.3,0.3),
        legend.title = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = NA),
  )+
  scale_y_reverse()+
  scale_x_reverse()+
  guides(colour = guide_legend(reverse = FALSE))

Figure_6<-ggarrange(TAXA_GK, AMR_GK, PA_erd_PCA,
                    font.label = list(size = 5),
                    common.legend = FALSE, nrow = 1)

pdf(file = "./R_output_files/Figures/Figure_6.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 4)
Figure_6
dev.off()


#Additional Matching Analysis####
Master$FreqEx<-ifelse(Master$ExacerbatorState %in% c("Exacerbator","NonEx"),"NFE","FE")
Master$FreqEx <- factor(Master$FreqEx, levels = c("NFE","FE")) #set ranking for order of bars

ALL<-ggplot(data=AMRFam,aes(x=Country, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('BC1', 'BC2'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=5), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRFam$Continent, scales="free_x", space = "free_x")+
  #facet_wrap(AMRFamMT$Country~AMRFamMT$FreqEx, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

TotalMatched<-ggplot(data=AMRFamMT,aes(x=Country, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('BC1', 'BC2'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=5), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRFamMT$Continent, scales="free_x", space = "free_x")+
  #facet_wrap(AMRFamMT$Country~AMRFamMT$FreqEx, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )
Paired<-ggplot(data=AMRFamPR,aes(x=Country, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('BC1', 'BC2'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=5), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRFamPR$Continent, scales="free_x", space = "free_x")+
  #facet_wrap(AMRFamMT$Country~AMRFamMT$FreqEx, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )


Trio<-ggplot(data=AMRFamTR,aes(x=Trio, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('BC1', 'BC2'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=5), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  #facet_grid(~AMRFamTR$Country, scales="free_x", space = "free_x")+
  #facet_wrap(AMRFamMT$Country~AMRFamMT$FreqEx, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

UnMatched<-ggplot(data=AMRFamUM,aes(x=Country, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('BC1', 'BC2'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=5), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRFamUM$Continent, scales="free_x", space = "free_x")+
  #facet_wrap(AMRFamMT$Country~AMRFamMT$FreqEx, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

GeoAet<-ggplot(data=AMRFamMT,aes(x=SC_AMR_alt, y=RPKM, fill=Resistome))+
  geom_bar(aes(), stat="identity", position = "fill") +
  scale_fill_manual(values = c("#026EB8","#06A955","#5D2E83","#2A2A73","#fc8403","#EBA5F3","#fc5017","#5CA5DB","#db6960","#a3d9d2","#B60004","#91CE59","#97809e","#C6DFA6","#FF9300","#FFBC06","#3B3B3B", "#026EB8","#06A955","#ffcccc","#2A2A73"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_discrete(labels = c('BC1', 'BC2'))+
  theme(legend.position="bottom",
        #axis.text=element_blank(),
        axis.title=element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))+ 
  guides(fill=guide_legend(nrow=5), size = .1)+
  xlab("")+
  ylab("Relative abundance (%)")+
  facet_grid(~AMRFamMT$Country, scales="free_x", space = "free_x")+
  #facet_wrap(AMRFamMT$Country~AMRFamMT$FreqEx, scales="free_x", nrow = 1)+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_text(size = 12)
  )

MatchFig<-ggarrange(TotalMatched, Paired, Trio, UnMatched,
                    font.label = list(size = 5),
                    common.legend = TRUE, legend = "bottom" ) 

chisq.test(Master$Aetiology_short, Master$Country)
table(Master$SC_AMR_alt)






