Fig_5<-ggarrange(p, p2, nrow=1)
pdf(file = "../Data/R_output_files/Fig_5DE.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 8)
Fig_5
dev.off()

pdf(file = "../Data/R_output_files/Fig_6new.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 6)
S7
dev.off()


pdf(file = "../Data/R_output_files/Fig_6new2.pdf",   # The directory you want to save the file in
              width = 6, # The width of the plot in inches
              height = 5)
PA_erd_PCA
dev.off()

Fig6new3<-ggarrange(TAXA_GK, AMR_GK, nrow=1)
pdf(file = "../Data/R_output_files/Fig6new3.pdf",   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 4)
Fig6new3
dev.off()
