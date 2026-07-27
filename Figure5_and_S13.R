library(tidyr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gghalves)
library(dplyr)
library(viridis)
library(ggnewscale)
library(tidyverse)
library(directlabels)
options(scipen=999) #non-scientific notation
ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12),
  axis.line = element_line(colour = "grey"),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10)
)

Dir="~/Paper/ModelSexChrom/V3/CleanDataset/"
ID="1750881723915" #Seed Id of the simulation to plot (figure 5) #Used to easily plot multiple simulation
#ID="1953648752606" #Seed Id of the simulation to plot (uncomment to create figure S24) 

Precision=10000 #Plot precision
Freq=intersect(intersect(list.files(Dir, pattern = ID, full.names = T),  #file containing inversion frequency
                         list.files(Dir, pattern = "Freq", full.names = T)),
               list.files(Dir, pattern = "parsed", full.names = T))
Recom=intersect(intersect(list.files(Dir, pattern = ID, full.names = T), #file containing recombination rate
                          list.files(Dir, pattern = "Recom", full.names = T)),
                list.files(Dir, pattern = "parsed", full.names = T)) 

Data=read.table(Recom, stringsAsFactors = F, header = F) #Read the recombination file (panel B)
Pos=seq(1:(length(colnames(Data))-1))
colnames(Data)=c("Generation", Pos) #Change col names with position in Mb
Data$"101"=Data$"101"-0.5 #The bin 101 contain the position that split the two chromosome in SLiM. This position recombine with a probability 0.5 each generation (as in nature). So remove these recombination event from the recombination rate
DataLong=gather(Data, Position, Value, "1":"200") #Pivot the plot
DataLong$Position=as.numeric(as.character(DataLong$Position)) #Define position as numeric for plotting
DataLong$Chrom="Chromosome 1 (X/Y)" #Change labels
DataLong[DataLong$Position>100,]$Chrom="Chromosome 2"
DataLong$NormVal=DataLong$Value/mean(DataLong[DataLong$Generation<10000,]$Value) #Define recombination rate relative to the recombination rate during the burn-in phase
DataLong[DataLong$Chrom=="Chromosome 2",]$Position=DataLong[DataLong$Chrom=="Chromosome 2",]$Position - 100 #Change the position of the autosome (from 101-200 to 1-100)
DataSub=subset(DataLong, (DataLong$Chrom=="Chromosome 1 (X/Y)" & DataLong$Generation < 101000)) #Only sex chrom data (ugly, just used for plotting the sex locus position)
base=ggplot(DataLong[DataLong$Generation < 101000,]) #Limit plot length
PlotRecomb=base+geom_tile(aes(x=Position, y=Generation, fill=NormVal),colour="white",size=0.02)+
  scale_y_reverse(expand=c(0.005,0),breaks=pretty(DataLong[DataLong$Generation < 101000,]$Generation, 10))+
  geom_vline(data=DataSub, aes(xintercept = 50), color="red2")+
  scale_color_manual("", values=c("red"))+
  scale_x_continuous(expand=c(0.01,0.01))+
  scale_fill_viridis("Relative \n recombination rate", option = "inferno", direction = -1, limits=c(0,max(DataLong$NormVal)))+
  facet_wrap(.~fct_relevel(as.factor(Chrom), "Chromosome 1 (X/Y)", "Chromosome 2"), scale="free_x")+
  xlab("Position (Mb)")+
  theme(panel.border = element_blank(),
        plot.margin = margin(0,0,0,0, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.title = element_text(face="bold"),
        axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=10),
        strip.text = element_blank())


DataInv=read.table(Freq, stringsAsFactors = F, header = F)  #Read the inversion frequency file (panel A)
colnames(DataInv)=c("Generation", "InvStart", "InvEnd", "InitMut", "p__", "pI_", "pII", "Yfreq", "Xfreq","MeanMutInv", "MeanMutNoInv")
DataInv$InvStart=DataInv$InvStart/1000000 #Define the inversion position on the new scale (in Mb instead of in bp)
DataInv$InvEnd=DataInv$InvEnd/1000000
DataInv$InvFreq=(DataInv$Xfreq + DataInv$Yfreq)/2000 #Overall inversion frequency (note that the Yfreq and Xfreq column contain the number of chromosome with inversion )
DataInv$Yfreq=DataInv$Yfreq/500 #Define inversion frequency on the Y and X chromosome (note their is n=2000 haploid genome, so 500 Y chromosomes)
DataInv$Xfreq=DataInv$Xfreq/1500
DataInv$Chrom="Chromosome 1 (X/Y)" #Define label
DataInv[DataInv$InvStart>100,]$Chrom="Chromosome 2"
DataInv[DataInv$Chrom=="Chromosome 2",]$InvStart=DataInv[DataInv$Chrom=="Chromosome 2",]$InvStart - 100 #Change the position of the autosome (from 101-200 to 1-100)
DataInv[DataInv$Chrom=="Chromosome 2",]$InvEnd=DataInv[DataInv$Chrom=="Chromosome 2",]$InvEnd - 100

Col=scales::viridis_pal(begin=0.00, end=0.80, option="cividis")(3)
DataInvSub=subset(DataInv, (DataInv$Chrom=="Chromosome 1 (X/Y)" & DataInv$Generation %% Precision == 0 & DataInv$Generation <101000))
base2=ggplot(DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000), ])
PlotFreq=base2+
  geom_rect(data=DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000 & DataInv$Chrom=="Chromosome 1 (X/Y)"), ],
            aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=Yfreq, fill="Y-linked"),color="black", alpha=0.6, size=0.2)+
  geom_rect(data=DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000 & DataInv$Chrom=="Chromosome 1 (X/Y)"), ], 
            aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=Xfreq, fill="X-linked"),color="black", alpha=0.6, size=0.2)+
  geom_rect(data=DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000 & DataInv$Chrom=="Chromosome 2"), ],
            aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=InvFreq, fill="Autosomal"),color="black", alpha=0.6, size=0.2)+
  scale_fill_manual("Inversions", values = Col)+
  ylab("Inversion Frequency")+
  scale_color_manual("", values=c("red"))+
  geom_vline(data=DataInvSub, aes(xintercept = 50), color="red2")+
  scale_x_continuous(expand=c(0.01,0.01), limits = c(0,100))+
  scale_y_continuous(expand=c(0,0),breaks = c(0.0,0.5,1.0), position = "right")+
  facet_grid(Generation~fct_relevel(as.factor(Chrom), "Chromosome 1 (X/Y)", "Chromosome 2"), scale="free_x", drop=F, switch="y")+
  theme(
    strip.text.y.left = element_text(angle = 0, size=10),
    strip.text.x = element_text(face="bold", size=14),
    panel.border = element_blank(),  
    plot.margin = margin(0,0,0,0,"pt"),
    panel.grid.major = element_line(color="grey", size=0.2),
    panel.grid.minor = element_line(color="grey", size=0.2),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(face="bold", size=12),
    legend.title = element_text(face="bold", size=12),
    axis.title.y = element_text(face="bold", size=14),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(size=8))

plots <- align_plots(PlotFreq, PlotRecomb, align = 'vh', axis = 'lrtb')

MergedPlot=plot_grid(plots[[1]],plots[[2]],  ncol=1, labels = c('a', 'b'))
Sub=sub(".*//", "",sub("_InvFreq.*", "", Freq))
save_plot(paste("~/Paper/ModelSexChrom/V3/Plot/Fig4.png"), MergedPlot, ncol = 2, nrow=2)
save_plot(paste("~/Paper/ModelSexChrom/V3/Plot/Fig4.pdf"), MergedPlot, ncol = 2, nrow=2)