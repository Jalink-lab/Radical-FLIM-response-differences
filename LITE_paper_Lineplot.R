#R workspace to plot lineplot  out of the .csv results from the Fiji macro
#SM & KJ June 2023 revised version for making lineplots

#Before running the code, these libraries need to be installed----------------------------------------------------------------------------------------------------------------------------------------
library("readxl")
library("ggplot2")
library("stringr")

#Loading and reading the file-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())                                                                                         #clear old environment
myFile <- file.choose()                                                                                 #choosing file of interest
df  <- read_xlsx(myFile)                                                                                #making it into a dataframe

df$Time_points <- factor(df$Time_points, levels = c("0s","1s","10s", "30s", "1min", "2min"))            #to plot the lineplot in this order

#Starting plotting here-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot<-ggplot(data=df, aes(x=Time_points, y=Mean_intensity, group=1))+
  
  geom_line(color="#009E73", size=1.2)+
  geom_point(color="darkgreen", size=3) +
  
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
  ylab("Mean intensity") + xlab("Time")


plot

line_plot <- plot
line_plot <- ggplot() +
  draw_plot(combi_plot, width = 1, height = 1)

line_plot

#Saving the plot--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_name <- "ROS_dye_intensity_lineplot"
plot_save_dir <- "D:/Local Surfdrive_SM/LITE_figures/Output_PDFs_20Traces"
plot_name <- paste(plot_name,".pdf")
ggsave(path = plot_save_dir, filename = plot_name, height=20, width=30, units=c("cm"), dpi=600)