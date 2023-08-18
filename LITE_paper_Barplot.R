#R workspace to plot barplot out of  the .csv results from the Fiji macro
#SM & KJ June 2023 revised version for making barplots

#Before running the code, these libraries need to be installed---------------------------------------------------------------------------------------------------------------------------------------
library("readxl")
library("ggplot2")
library("stringr")

#Loading and reading the file------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())                                                                                               #clear old environment
myFile <- file.choose()                                                                                       #choosing file of interest
df  <- read_xlsx(myFile)                                                                                      #making it into a dataframe

df$Sample <- factor(df$Sample, levels = c("Unstimulated","40nM IsoP","40nM IsoP_Blue light 430nm(2 mins)"))   #to plot the barplot in this order

#Starting plotting here------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot<-ggplot(data=df, aes(x=Sample, y=Concentration))+
            
   geom_bar(aes(color = Replicate, fill = Replicate),
              stat = "identity", position = position_dodge(0.8),width = 0.7) +
              scale_color_manual(values = c("#56B4E9", "#009E73"))+                                           
              scale_fill_manual(values = c("#56B4E9", "#009E73")) +

theme_classic(base_size = 16) + theme(panel.background = element_blank(),                                     
                                      panel.grid.major = element_blank(),                                     
                                      panel.grid.minor = element_blank(),legend.position =  "none",
                                      axis.line = element_line(colour = "black", size = 1.6)) +
  theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size = 20, color = "black") ,
        axis.text.y = element_text(size = 20, color = "black"),
        legend.title=element_text(size=20), 
        legend.text =element_text(size=20)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) + 
  ylab("cAMP concentration(nM)") + xlab("")


plot

bar_plot <- plot
bar_plot <- ggplot() +
  draw_plot(combi_plot, width = 1, height = 1)

bar_plot

#Saving the plot-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_name <- "cAMP_ELISA_Barplot_Cos7WT"
plot_save_dir <- "D:/Local Surfdrive_SM/LITE_figures/Output_PDFs_20Traces"
plot_name <- paste(plot_name,".pdf")
ggsave(path = plot_save_dir, filename = plot_name, height=20, width=30, units=c("cm"), dpi=600)



