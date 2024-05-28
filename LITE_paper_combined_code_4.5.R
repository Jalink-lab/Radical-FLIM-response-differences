# R workspace to plot individual flim traces out of the .csv results from the Fiji macro
# SM & KJ April 2023 revised version for batch processing an entire folder of locFiles.
# the script plots a number of cells (nrOfTraces) randomly picked from all cells that obey to
# baseline = within 'baseLineRange' ns from the median value, and excludes traces that have many NaN.
# it also plots measurements from each of these cells taken at timepoints between the two gray timelines,
# along with some statistics of those.

#Before running the code, these libraries need to be installed---------------------------------------------------------------------------------------------------------------------------------------
library("tidyr")
library("ggplot2")
library("cowplot")
library("patchwork")
library("RColorBrewer")
library("pals")
library("dplyr")

#Definition of main function-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SM_func1<- function(locFile, locNrOfTraces=20, locBaseLineCutOff=10, locBaseLineRange=0.1, locTraceEndTime=-1, locScatterBoxPlotTime=-1, locScatterBoxPlotRunLength=3) { 
  #S.M. March 2023. this version contains all code to read a single csv file and make a combiplot of it. 
  #use it in a loop to plot all files. 
  df  <- read.csv(locFile, sep = ',')                 #next, first curate df in case any NAN
  df <- df[ , colSums(is.na(df))==0]                  #removes columns which have any NaNs
  if (any(is.na(df))){    
    print(" WARNING: CSV file contains NaN values. Please curate that file outside this routine!! Processing skipped.")
  }
  nCells = ncol(df)-1                                 #verify nr of cells and nr of timepoints
  if(locNrOfTraces < 0 | locNrOfTraces > nCells){     #plot locNrOfTraces cells, all traces if -1, or maximally nCells-1
    locNrOfTraces = nCells
  }
  nSamples = nrow(df)
  timePerFrame = (df[nSamples,1]-df[1,1])/(nSamples-1) 
  
  colnames(df)[1]<-"Time"                             #Changing column name to 'Time' - will be useful later for tidying the data
  colnames(df)[2]<-"Y0"                               #Change name of the first cell 
  
  if (locScatterBoxPlotTime < 0){
    locScatterBoxPlotTime=df[nSamples,1]             
  }
  print(locScatterBoxPlotTime)                        #print scatterboxplottime in console
  
  mean_baseline_value <- median(as.matrix(kj<-df[1:locBaseLineCutOff,2:nCells]))  #mean baseline values for the first 'locBaseLineCutOff' time frames has been extracted
  pos_range <- mean_baseline_value + locBaseLineRange
  neg_range <- mean_baseline_value - locBaseLineRange   
  
  #next block: Clearing out the data with poor baseline - here baseline values +/-0.1 of mean baseline value is kept, rest is discarded-----
  columns_to_check <- names(df)[-1]                  # get the names of all columns except the first one(which is Time)
  df_copy <- df
  for(col in columns_to_check) {                     # go over all columns and check the values in the "locBaseLineCutOff" value)
    test_TF <- between(abs(df[1:locBaseLineCutOff, col]), neg_range, pos_range) #testing if the condition will be True or False
    test_all <- all(test_TF)                         # test if within range
    if (test_all == FALSE)                           # reject columns that are out of range
    {
      df<- select(df, -c(col))
    }
    rm(test_all)                                     #remove these variables to prevent overlaps
    rm(test_TF)                                      #remove these variables to prevent overlaps
  }
  locNrOfTraces<-min(locNrOfTraces, (ncol(df)-1))    #maximum nr of traces is all of them
  nLabel=paste(ncol(df),"cells")
  
  # next block: choose locNrOfTraces random column names-----
  columns_to_include <- sample(names(df)[-1], locNrOfTraces)    # selecting the columns of 'n' random columns 
  new_df <- df[, columns_to_include]                            # make new df with the selected columns
  new_df <- cbind(Time = df$Time, new_df)                       # add the "Time" column to the new dataframe 
  new_df_tidy <- gather(new_df, Cell, Lifetime, -Time)          # tidying the data
  end_point <- df_copy[nrow(df_copy),1]                         # calculating automatically where the traces end
  
  #Next block: plot those random single cells in colors-----
  locTraceEndTime<-max(locTraceEndTime, df[nSamples,1])         #if traceEndTime ==-1 then use last sample
  my_colors <- brewer.pal(12, "Paired")                         #Choosing a color palette 
  my_colors <- colorRampPalette(my_colors)(locNrOfTraces)       # Add more colors to this palette :
  
  #Starting the single cell traces plotting here:
  plot <- ggplot() +
    geom_line(data =  new_df_tidy, mapping = aes(x=Time, y=Lifetime, color = Cell), size= 1, alpha=0.8)+
    scale_color_manual(values = my_colors) +
    coord_cartesian(ylim = c(2.25, 3.4), xlim = c(0,end_point))+
    
    theme_classic(base_size = 16) + theme(panel.background = element_blank(),                                
                                          panel.grid.major = element_blank(),                                
                                          panel.grid.minor = element_blank(),legend.position =  "none",
                                          axis.line = element_line(colour = "black", size = 1.6)) +
    theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"),
          axis.title = element_text(size=20,face="bold"),
          axis.text.x = element_text(size = 20, color = "black") ,
          axis.text.y = element_text(size = 20, color = "black"),
          legend.title=element_text(size=20), 
          legend.text =element_text(size=20)                                                        
    ) +
    geom_vline(xintercept = locScatterBoxPlotTime, linetype="dotted", color = "gray", size=1) +
    geom_vline(xintercept = (locScatterBoxPlotTime-timePerFrame*locScatterBoxPlotRunLength), linetype="dotted", color = "gray", size=1) +
    annotate("text", x=0.85*locTraceEndTime, y=2.26, label= nLabel, size =6) + 
    ylab("Lifetime (ns)") + xlab("Time (s)")
  
  #  print(plot)
  
  #Starting the boxplot plotting here:
  scatterBoxPlotSample=round(locScatterBoxPlotTime/timePerFrame)
  df_boxplot<-tail((head(df, scatterBoxPlotSample-locScatterBoxPlotRunLength)), locScatterBoxPlotRunLength) #else locScatterBoxPlotRunLength samples taken starting from the specified scatterBoxPLotTime
  
  df_boxplot <- data.frame(colMeans(df_boxplot[,-1]))
  colnames(df_boxplot)[1] <- "X"                                                                           #contains means of 'scatterBoxPLotRunLength' samples taken starting at 'scatterBoxPLotTime'
  
  boxplot = ggplot(df_boxplot, aes(x = "", y = X)) + 
    
    stat_boxplot(geom = "errorbar",
                 width = 0.3,               
                 alpha = 1,                 
                 size = 0.6) +              
    
    geom_boxplot(width = 0.6,                  
                 alpha = 1,                 
                 color = "black",           
                 outlier.shape = NA,         
                 lwd = 0.6, 
                 fatten = 0.6)          +   
    
    coord_cartesian(ylim = c(2.25, 3.4)) +
    
    geom_point(position = position_jitter(seed = 0.5, width = 0.15), size =0.2, color ="#FF7F00") +
    
    theme_void(base_size = 16) + theme(aspect.ratio =5.5)
  
  boxplot 
 
  
  #Plot both plots together
  combi_plot <- plot + boxplot
  combi_plot <- ggplot() +
    draw_plot(combi_plot, width = 1, height = 1)
  
  combi_plot
  
  return(combi_plot)
  
}
#end def main function


#set up data path------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wd<-"D:\\Local Surfdrive_SM\\LITE_figures"     #path to the folder that contains both csv folder and output folder(hard-coded; replace with your own path)
setwd(wd)
dir<-paste(wd,"\\INPUT_CSVS_FINAL",sep="")     #sep has to be "" because default is " "
plot_save_dir <- paste(wd,"\\Output_Temp",sep="")


#FILE LOOP BLOCK for running all files. get a list of all files with csv in the name in your directory-----------------------------------------------------------------------------------------------
files<-list.files(path=dir, pattern='.csv', full.names = TRUE)
Filenames_changed <- gsub("/", "\\", files[1:length(files)], fixed=TRUE)

for (k in 1:length(Filenames_changed)){
  print(Filenames_changed[k])
  plot_name <- substr(Filenames_changed[k], (nchar(dir)+2), nchar(Filenames_changed[k])-4)
  scatterBoxPlotTime=-1
  if(substr(plot_name,1,1)=="T"){
    scatterBoxPlotTime=as.numeric(substr(plot_name, 2, 5) )        # grab the time at which to scatterplot from the filename
  }  
  print(plot_name)
  SM_func1(Filenames_changed[k], nrOfTraces, baseLineCutOff, baseLineRange, traceEndTime, scatterBoxPlotTime, scatterBoxPlotRunLength)

  plot_name <- paste(plot_name,".pdf")
  ggsave(path = plot_save_dir, filename = plot_name, height=20, width=30, units=c("cm"), dpi=600)
  rm(plot_name)
  scatterBoxPlotTime=-1 # in case you want to make single plots after this
}
#END OF FILE LOOP BLOCK


#TEMP BLOCKS FOR CODE ERROR TESTINGS-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#TEMP BLOCK: CODE OVERRIDE PARAMS to make testing easier
nrOfTraces=-1                             #-1 signals: plot all
baseLineCutOff=10                         #default = 10. Nr of samples to average for detecting baseline
baseLineRange=0.1                        #0.1 ns default
traceEndTime=-1                          #-1 signals: plot end of trace
scatterBoxPlotTime=-1                   #option to set timepoint (-1 = end) at which the END of the stretch of values for scatterplot is taken
scatterBoxPlotRunLength=5              #default = 3

#---alternative code snippets. select necessary lines to run a single file
setwd(dir)
#File <- file.choose()

#SM_func1(File, nrOfTraces, baseLineCutOff, baseLineRange, traceEndTime, scatterBoxPlotTime, scatterBoxPlotRunLength) 
#END OF TEMP OVERRIDE BLOCK. remove later


#TEMP block for error seeking
locFile=File  #this makes the local parameters global so they show up in Environment.
LocNrOfTraces=nrOfTraces
locBaseLineCutOff=baseLineCutOff
locBaseLineRange=baseLineRange
locTraceEndTime=traceEndTime
locScatterBoxPlotTime=scatterBoxPlotTime=
locScatterBoxPlotRunLength=scatterBoxPlotRunLength  
#end error seeking block


