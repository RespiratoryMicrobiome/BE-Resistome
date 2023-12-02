library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(vegan)

###PERMANOVA CHECKS
AMRLT_gene <- MasterLT %>%
  as_tibble() %>%
  select(1:13,80:230) 
AMR_cols<-colnames(AMRLT_gene[14:164])
AMRLT_gene <- AMRLT_gene %>%
  gather(AMR, RPKM, AMR_cols, -SampleSeqNo, -SputumSampleNo,  -TypeSamples, -TypeSamplesA,-TypeSamplesB,-Exacerbations,	-FEV1, -BSI,	-Severity, -WksToNxtEx,	-TmToNxtEx,	-Antibiotic,	-Antibiotic_class)
AMRLT_gene$TmToNxtEx <- factor(AMRLT_gene$TmToNxtEx , levels = c("MoreThan12w","LessThan12w"))
AMRLT_gene$Exacerbations <- factor(AMRLT_gene$Exacerbations , levels = c("NFE","FE"))
relapse.labs <- c(
  `LessThan12w` = "<12 w",
  `MoreThan12w` = ">12 w")
AMRLT_gene$FEV170<-ifelse(AMRLT_gene$FEV1 >70, ">70", "<70")
AMRLTctrols<-subset(AMRLT_gene, is.na(TypeSamplesB))


AMRLT_diversity <- MasterLT[which(MasterLT$TypeSamplesA !="NA"),] %>%
  as_tibble() %>%
  select(1,80:230) 

main_dataLT <- AMRLT_diversity
AMRLT_diversity<-as.matrix(AMRLT_diversity)
rownames(AMRLT_diversity) <- AMRLT_diversity[,1]
AMRLT_diversity = as.data.frame(subset(AMRLT_diversity, select = -c(SampleSeqNo) ))
AMRLT_diversity[] <- lapply(AMRLT_diversity, as.numeric)
isZero <- base::rowSums(AMRLT_diversity) == 0
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
E<-ggplot(gg) +
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


#PERMANOVA - timepoint
adonis2(AMRLT_diversity~TypeSamplesB , data = LTDiversityViz, method = "bray",permutations=999, strata = LTDiversityViz$SputumSampleNo)

#Exclude samples with 'EX' in TypeSamplesB so as to compare pre post ACUTE therapy
filtered_data <- LTDiversityViz[LTDiversityViz$TypeSamplesB != 'EX', ]

filtered_data <- LTDiversityViz[LTDiversityViz$TypeSamplesA == "P1" | LTDiversityViz$TypeSamplesA == "BSL",]

# Assuming AMRLT_diversity is your dataset with values
# Match the rows based on SampleSeqNo
amr_diversity <- AMRLT_diversity[rownames(AMRLT_diversity) %in% filtered_data$SampleSeqNo, ]

# Run PERMANOVA with strata information
result <- adonis2(amr_diversity ~ TypeSamplesB, data = filtered_data, method = "bray", 
                  permutations = 999, strata = filtered_data$SputumSampleNo)


# Print the result
print(result)







