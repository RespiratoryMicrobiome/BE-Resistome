library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)

setwd("c:/Users/mmaca/Code/git/BE-Resistome/Analysis/")
#Load data
qPCR <-read.csv("../Data/R_input_files/qPCR_data.csv") %>%
  as_tibble()

long_qPCR <- gather(qPCR, key = "Copies", value = "Value", PA_qPCR, HI_qPCR, X16S_copies)

long_qPCR$Value_log <- log10(long_qPCR$Value + 0.1)

long_qPCR$Copies<- factor(long_qPCR$Copies, levels = c("X16S_copies", "HI_qPCR", "PA_qPCR"))

qPCR_plot <-ggplot(long_qPCR, aes(x = Copies, y = Value_log, color = Copies)) +
  geom_boxplot() + 
  scale_color_manual(values = c("HI_qPCR" = "#BEAED4", "PA_qPCR" = "#7FC97F","X16S_copies" = "#555555" ), labels = c("HI_qPCR" = "HI", "PA_qPCR" = "PA", "X16S_copies" = "16S (total)"))+
  geom_jitter(width = 0.1) + # You can change the type of plot based on your preference
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8,10))+
  facet_wrap(~ RT_group, scales = "free_y", labeller = labeller(RT_group = c("1" = "RT1", "2" = "RT2")))+
  labs(x = "", y = expression(paste("log10 (16S rRNA gene copies /", mu, "l)"))) +
  #scale_y_log10() +  # Add this line to set the y-axis to log scale
  theme_minimal()+
  theme(axis.line = element_line(size = 0.5, colour = "black"))+
  scale_x_discrete(labels = c("16S (total)", "HI", "PA"))+
  guides(
    color = guide_legend(
      title = NULL)) 

# qPCR_plot + stat_compare_means(
#   aes(group = RT_group),
#   comparisons = list(c("HI_qPCR", "PA_qPCR")),
#   method = "wilcox.test",
#   label = "p.signif",  # Use "p.signif" to display significance stars
#   test.adj = "bonferroni",
#   position = "identity",
#   vjust = -0.5
# )


# Calculate predominance based on which value is higher
df <- qPCR[!grepl("GREEK", qPCR$Sample.Name, ignore.case = TRUE), ] %>%
  mutate(Predominance = ifelse(PA_qPCR > HI_qPCR, "PA_qPCR", "HI_qPCR"))

table(df$RT_group, df$Predominance)

# Calculate predominance based on which value is higher
df <- df %>%
  mutate(Predominance = ifelse(PA_qPCR > HI_qPCR, "PsA-dom", "Hi-dom"))

filtered_df <- df %>%
  filter(RT_group %in% c("1", "2"))

# Calculate proportions within each group
proportions_df <- filtered_df %>%
  group_by(RT_group, Predominance) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(RT_group) %>%
  mutate(Proportion = Count / sum(Count))

qPCR_prop<- ggplot(proportions_df, aes(x = as.factor(RT_group), y = Proportion, fill = Predominance)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Proportion", fill = "Predominance") +
  scale_fill_manual(values = c("PsA-dom" = "#7FC97F", "Hi-dom" = "#BEAED4")) +
  scale_x_discrete(position = "top") +
  theme_minimal()


combined_plot <- ggarrange(qPCR_plot, qPCR_prop, ncol = 2, #labels = "AUTO",
                           common.legend = FALSE,
                           legend = "bottom", widths = c(1.5, 1))

pdf(file = "../Data/R_output_files/Fig_E9.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6)
combined_plot
dev.off()


#Load data
qPCR <-read.csv("../Data/R_input_files/qPCR_data_R.csv") %>%
  as_tibble()

qPCR_2 <- qPCR[grepl("GREEK", qPCR$Sample.Name, ignore.case = TRUE), ]

long_qPCR_2 <- gather(qPCR_2, key = "Copies", value = "Value", PA_qPCR, HI_qPCR, X16S_copies)

long_qPCR_2$Value_log <- log10(long_qPCR_2$Value + 1)

long_qPCR_2$Copies<- factor(long_qPCR_2$Copies, levels = c("X16S_copies", "HI_qPCR", "PA_qPCR"))
long_qPCR_2$Time.point<- factor(long_qPCR_2$Time.point, levels = c("Pre", "Post"))

qPCR_plot_PA_erd <-ggplot(long_qPCR_2[which(long_qPCR_2$Copies == "PA_qPCR"),], aes(x = Time.point, y = Value_log)) +
  geom_boxplot(outlier.shape = NA, fill = "#7FC97F")+ 
  geom_line(aes(group = Patient))+ 
  geom_point() + # You can change the type of plot based on your preference
  #facet_wrap(~ Copies, scales = "free_y", labeller = labeller(RT_group = c("1" = "RT1", "2" = "RT2")))+
  labs(x = "", y = expression(paste("log10 (16S P.aeruginosa rRNA gene copies /", mu, "l)"))) +  scale_x_discrete(labels = c("Pre", "Post"))+
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4,5,6))+
  #scale_y_log10() +  # Add this line to set the y-axis to log scale
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1, linetype="solid"),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = NA),
    axis.line = element_line(size = 0.5, colour = "black"))

wilcox_result<-
wilcox.test(long_qPCR_2[which(long_qPCR_2$Copies== "PA_qPCR"),]$Value~long_qPCR_2[which(long_qPCR_2$Copies== "PA_qPCR"),]$Time.point, paired = TRUE, alternative = "two.sided")

pdf(file = "../Data/R_output_files/Fig_E7.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6)
combined_plot
dev.off()

pdf(file = "../Data/R_output_files/Figure_E11.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
qPCR_plot_PA_erd
dev.off()









