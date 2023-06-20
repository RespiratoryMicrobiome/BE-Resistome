#SETUP####
#Load packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, ggplot2, tidyverse, lubridate, data.table, tsibble, wesanderson, reticulate, SNFtool, funrar, vegan, dunn.test, ggpubr, Hmisc, RColorBrewer, pheatmap, phyloseq, ggpubr, colorspace, RColorBrewer, grid)
source("R_input_files/function_snf.R")
use_python("/usr/local/lib/python3.7/") #something up with reticulate
source_python("R_input_files/sil.py")
#py_install("scikit-learn")

#Ivan set parameters####
project_name <- "CAMEB2"
wkdir <- file.path("/Users/michealmacaogain/Dropbox/Research/Publications/44. CAMEB_2/Ivan_AMR_analysis_2020/", project_name)
setwd(wkdir)
Omic <- "Microbiome-AMR"
groupvar <- "SC_AMR_alt"
#groupvar <- "Continent"
my.labels <- c("H", "B")
#my.labels <- c("0", "1", "2")

# Create output directory
outdir <- file.path(wkdir, "analysis", "HEATMAP")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE) 

#Species####
taxLevel <- "Gene"

##Read ps object####
ps <- readRDS(file.path(wkdir, "analysis", "AMR", paste0("ps_", taxLevel, ".RData")))

#.AMR data####
ab_data=make_relative(as.matrix(read.csv("../../Take_over/MASTER DATA/Spectral_clustering/AMR_incl_helathyctrl_Drug.csv", row.names = 1)))*100 #load functional data
ab_data[is.nan(ab_data)] <- 0
ab_data<-as.data.frame(ab_data)
#filter based on prevalence
#z=colSums(ab_data>0.1) #filter to reduced number of taxa 
#sel_col=row.names(as.data.frame(z[z>=(0.01*(nrow(ab_data)))])) #In 1% patients prevalent
#ab_data<-ab_data[sel_col]
#remove(sel_col,z)
#ab_data<-ab_data[rowSums(ab_data[, -1])>0, ] #drop no_res samples
#ab_data<-ab_data[row.names(ab_data) != "TBS672", , drop = FALSE] #TBS672 contains a single PatA gene creating clustering artifact

# ab_data<-as.data.frame(make_relative(as.matrix(ab_data)))
# write.csv(ab_data, "R_output_files/AMRFiltered.csv")
# ab_data

#Master data####
Master <-read.csv("../../Take_over/MASTER DATA/R_input_files/MetaDataWsnfResults_3.5temp_TJ.csv") %>%
  #Master <-read.csv("R_input_files/MetaDataWsnfResults_3.5temp_TJ.csv") %>% # v 3.3 includes raw read data, 3.4 includes info on paring and "trios"
  as_tibble()
Master$PseuDOM<-ifelse(Master$DomBactOth_short == "Pseudomonas", "Pseudomonas", "other")
Master$HaemDOM<-ifelse(Master$DomBactOth_short == "Haemophilus", "Haemophilus", "other")
Master$FEVfactor<-cut(Master$FEV1, breaks=c(0, 30, 50, 70, Inf))
Master$AetiologyN3<-ifelse(Master$Aetiology_short == "postTB", "postInfect", Master$Aetiology_short)
Master$SampleID <- factor(Master$SampleID, levels = Master$SampleID[order(Master$SC_V)]) #set ranking for order of bars
#subset on matched patients
MasterMT<-subset(Master, Matching == "Matched")
MasterPR<-subset(Master, Paired == "Pair")
MasterTR<-subset(Master, Trio != "NonTrio")
MasterSick<-subset(Master, SC_AMR_alt != "NA")
#[which,(MasterSick$SC_AMR_alt != "NA"),]

##Sample ordering####
ab_data$SC_AMR_alt<-ifelse(is.na(Master$SC_AMR_alt) == TRUE,  "H", Master$SC_AMR_alt)
#ab_data$SC_AMR_alt<-ifelse(is.na(Master$SC_AMR_alt) == TRUE,  "H", "B")
#ab_data$SC_AMR_alt <- factor(ab_data$SC_AMR_alt, levels = c("H","B")) #set ranking for order of bars

groupvar = "SC_AMR_alt"

sorted_group <- sort(as.vector(data.frame(ab_data))[ , groupvar], index.return= TRUE)

starts <- which(!duplicated(sorted_group$x))
sample_order <- c()
for (i in 1:length(starts)) {
  if (i!=length(starts)) {
    sample_order <- c(sample_order, sorted_group$ix[starts[i]:(starts[i+1]-1)])
  } else {
    sample_order <- c(sample_order, sorted_group$ix[starts[i]:length(sorted_group$ix)])
  }
}

#.AMR data reload####
ab_data_otu=make_relative(as.matrix(read.csv("../../Take_over/MASTER DATA/Spectral_clustering/AMR_incl_helathyctrl_Drug.csv", row.names = 1)))*100 #load functional data
ab_data_otu[is.nan(ab_data_otu)] <- 0
ab_data_otu<-as.data.frame(ab_data_otu)

##Prepare heatmap components####

lefsetaxa<-(ab_data_otu)
df <- t(ab_data_otu[sample_order, colnames(lefsetaxa)])
#rownames(df) <- c(tax_table(ps.prop)[rownames(df), taxLevel])
df <- log(df+1,2)

mymat <- as.matrix(data.frame(df))
mydf <- data.frame(row.names= as.vector(data.frame(Master)[sample_order, "SampleID"]), category= as.vector(data.frame(Master)[sample_order, groupvar]))

##Color and text setting####
#lightenParams <- seq(0, 0.9999, 1/55)
#col1 <- lighten("#152CAD", lightenParams)[length(lightenParams):1]
col1 <- sequential_hcl(12, h = 150, c = 80, l = c(35, 95), rev = TRUE, power = 2)[c(3, rep(4, 2), rep(5,3), rep(6,4), rep(7,5), rep(8,6), rep(9,7), rep(10,8), rep(11,9), rep(12, 10))]
newnames <- lapply(rownames(mymat), function(x) bquote(italic(.(x))))
my.labelsCol <- c("#E6B6FF", "#DF9A8F", "#1B9E77", "cyan"); names(my.labelsCol) <- my.labels
ann_colors = list(`category` = my.labelsCol)

##Plotting heatmap####
pheatmap(mymat, color= col1, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = mydf, annotation_colors = ann_colors,
         show_colnames = FALSE, labels_row = as.expression(newnames), cellheight = 8, cellwidth = 3,
         gaps_col = starts[2:length(starts)]-1, breaks= seq(0, 4.5, 0.1),
         annotation_names_row= FALSE, annotation_names_col= FALSE, annotation_legend= FALSE, fontsize= 7,
         file = file.path(outdir, paste0(groupvar, "_", taxLevel, "-heatmap.png")))

