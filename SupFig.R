## Sup Figures ##
library(tidyverse)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggraster)
library(grid)
library(gridExtra)
library(ggpubr)
library(Hmisc)
library(patchwork)
library(scales)
library(ggh4x)
library(ggrastr)
library(cowplot)
library(directlabels)
library(data.table)
### Figure S1 ###


Line=paste("time","h","s","u","n","r","Recomb","P",
           "FXN","FXI","FXNm","FXIm","FXNf","FXIf","FXm","FXf",
           "FYN","FYI","FY","Wm","Wf","D","q", sep=" ")
u=1e-08
n=2000000
File=paste0("~/Project/MutationShelteringV2/Revision1/Output/TimeSimul_h_s_u=",u, "n=", n, "_XYsyst_NoMutAccumul_Example.txt")
fwrite(list(Line), File) #Directly write the result in a file, for computing purpose
for (Recomb in c(0.0,0.5)) #Linked or unlinked
{
  for (h in c(0.1))
  {
    for (s in c(0.001,0.01,0.1))
    {
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      P=0.95 #Inversion with 5% less mutation than average
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXNm=1.00 #Frequency of X chromosomes in males with the non-inverted segment
      FXIm=0.00  #Frequency of X chromosomes in males with the inverted segment
      FXNf=1.00  #Frequency of X chromosomes in females with the non-inverted segment
      FXIf=0.00 #Frequency of X chromosomes in females with the inverted segment
      FYN=0.99 #Frequency of Y chromosome with the non-inverted segment
      FYI=0.01 #Frequency of Y chromosome with the non-inverted segment #Introduce the inversion in 1% of Y chromosomes
      FY=FYN+FYI #Frequency of Y chromosome (must equal 1, used for check)
      FXm=FXNm+FXIm  #Frequency of X chromosome in males(must equal 1, used for check)
      FXf=FXNf+FXIf #Frequency of X chromosome in females (must equal 1, used for check)
      FXI=(2/3)*FXIf + (1/3)*FXIm #Two third of the X chromosome are in females and one third in males
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII # Mean fitness of the males
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII  #mean fitness of the females
      D=FXI*FYN - FXN*FYI #Linkage disequilibrium
      time=0 #Initial time
      Line=paste(time,h,s,u,n,m,Recomb,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q,sep=" ") #Table row
      fwrite(list(Line), File, append = T)
      for (time in seq(2,10000,1)) #During 10000 generation
      {  #Recalculate the frequency of each inversion in each chromosome (see appendix for detail)
        FYI_t=(FYI*FXIf*WII +
                 FYI*FXNf*WNI*(1-Recomb) + 
                 FYN*FXIf*WNI*Recomb)/Wm
        FYN_t=(FYN*FXNf*WNN +
                 FYN*FXIf*WNI*(1-Recomb) + 
                 FYI*FXNf*WNI*Recomb)/Wm
        FXIm_t=(FXIf*FYI*WII +
                  FXIf*FYN*WNI*(1-Recomb) + 
                  FXNf*FYI*WNI*Recomb)/Wm
        FXNm_t=(FXNf*FYN*WNN +
                  FXNf*FYI*WNI*(1-Recomb) + 
                  FXIf*FYN*WNI*Recomb)/Wm
        FXIf_t=(FXIf*FXIm*WII +
                  (1/2)*(FXIm*FXNf + FXIf*FXNm)*WNI)/Wf
        FXNf_t=(FXNf*FXNm*WNN +
                  (1/2)*(FXIm*FXNf + FXIf*FXNm)*WNI)/Wf
        FYI=FYI_t #For the next generation, define the new inversion frequency
        FYN=FYN_t
        FXIm=FXIm_t
        FXNm=FXNm_t
        FXIf=FXIf_t
        FXNf=FXNf_t
        FXI=(2/3)*FXIf + (1/3)*FXIm
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Line=paste(time,h,s,u,n,m,Recomb,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q,sep=" ")
        fwrite(list(Line), File, append = T)
      }
    }
  }
}

Table=read.table(File, sep=" ", header=T, stringsAsFactors = F) #Read the dataset
Table$s=-Table$s
Table$Linkage="Autosomal" #Add labels
Table[Table$Recomb=="0",]$Linkage="Y-linked"
Table$FI=0.25*Table$FYI + 0.75*Table$FXI #Frequency of the inversion when not linked to the Y chromosome
Table[Table$Linkage=="Y-Linked",]$FI=Table[Table$Linkage=="Y-Linked",]$FYI #Frequency of the inversion when linked to the Y chromosome (frequency among the Y chromosome)

# Panel A #
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color palette
base=ggplot(Table)
PlotTrajDeter=base+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  geom_line(aes(x=time, y=FYI, linetype=as.factor(s), color=Linkage), size=1, alpha=1.0)+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_linetype_manual("s",values=c("solid", "dashed", "dotted"))+
  scale_color_manual("",values=Col, labels=c("Autosome", "Y"))+
  labs(title = "Deterministic trajectory of inversions")+
  ThemeSobr+
  theme(panel.border = element_blank(),  
        legend.key.width = unit(1.2,"cm"),
        legend.background = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(size=18, face="bold",hjust=0.5),
        strip.placement = "outside",
        axis.title = element_text(size=18),
        text = element_text(size=16),
        axis.line = element_line(colour = "grey"))+
  guides(linetype = guide_legend(title.position = "left", 
                                 label.position = "top",
                                 label.vjust = -2,
  ))+
  scale_y_continuous(breaks=c(0.0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.01,1.05))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))#+

PlotTrajDeter

# Panel B #
ResultHeatNonNeutral=read.table("~/Project/MutationShelteringV2/Revision1/Output/DeterministicSimul_sqrt.tsv", header=T, stringsAsFactors = F)
ResultHeatNonNeutralLong=ResultHeatNonNeutral %>% pivot_longer(cols=c("FracFixedY", "FracFixedAuto"), values_to = "FracFixed", names_to = "Chrom")

ResultHeatNonNeutralLong[ResultHeatNonNeutralLong$Chrom=="FracFixedY",]$Chrom="Y"
ResultHeatNonNeutralLong[ResultHeatNonNeutralLong$Chrom=="FracFixedAuto",]$Chrom="Autosome"

ResultHeatNonNeutralLong$s=-ResultHeatNonNeutralLong$s
base=ggplot(ResultHeatNonNeutralLong[ResultHeatNonNeutralLong$U %in% c(0.001, 0.01, 0.05),])
# Panel A
signed_sqrt_trans <- function() {
  trans_new(
    "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse   = function(x) sign(x) * (x^2),
    domain = c(-Inf, Inf)
  )
}

options(scipen = 999)
FracFixed=base+geom_tile(aes(x=h, y=s, fill=FracFixed))+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  # coord_cartesian(ylim=c(0,0.1), expand = F)+
  facet_grid(
    paste0("U[L]==", U) ~ Chrom,
    labeller = label_parsed
  )+
  scale_y_continuous(trans = signed_sqrt_trans(), breaks=c(-0.001, -0.01, -0.1))+
  scale_fill_viridis(option = "B", "")+
  ThemeSobr+
  theme(panel.border = element_blank(),
        legend.text = element_text(face="bold", size=14),
        strip.text.x = element_text(face="bold", size=16),
        strip.text.y = element_text( size=16),
        text = element_text(size=16),
        axis.title = element_text(size=18),
        legend.title = element_text(face="bold"),
        legend.background = element_blank(),
        plot.title=element_text(size=18, face="bold",hjust=0.5),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        #plot.margin = margin(29, 3, 3, 5, "pt"),
  )+ggtitle("Fraction of inversions expected to fix without drift \n and mutation accumulation")
FracFixed

PlotDeter=PlotTrajDeter / FracFixed + plot_annotation(tag_levels = list(c("a", "b"))) &
  theme(plot.tag = element_text(size = 16, face="bold"))
save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/PlotDeter.pdf"), PlotDeter, ncol=2, nrow=3)    

####################################################

### Figure S2 and S3 ###

options(scipen = 0)
Simul=read.table(paste("~/Project/MutationShelteringV2/Revision1/Output/IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat.24Col.txt",sep=""), stringsAsFactors = F) # Data for N=10000
colnames(Simul)=c("N", "u", "r", "h", "s", "sInv", "Gen", "StartInv", "EndInv", "Rep" ,"Chromosome",
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")
SimulCp=Simul
Simul$Gen=Simul$Gen - 6*Simul$N
Simul[Simul$Chromosome=="SexChrom",]$Freq=Simul[Simul$Chromosome=="SexChrom",]$Freq * 4 # Define Y inversion frequency
Simul$InvSize=Simul$EndInv - Simul$StartInv #Inversion size
Simul[Simul$Chromosome=="SexChrom",]$Chromosome="Ychromosome"

## Load the file containing the simulation results ##

# Figure S2 #

Simul1=Simul[(Simul$u==1e-08 & Simul$h==0.2 & Simul$s==-0.001 & Simul$InvSize==5000000 & Simul$N==1000 & Simul$Chromosome=="Ychromosome"),]
Simul2=Simul[(Simul$u==1e-08 & Simul$h==0.1 & Simul$s==-0.01 & Simul$InvSize==1000000 & Simul$N==1000 & Simul$Chromosome=="Ychromosome"),]
Simul3=Simul[(Simul$u==1e-08 & Simul$h==0.01 & Simul$s==-0.1 & Simul$InvSize==5000000 & Simul$N==10000 & Simul$Chromosome=="Ychromosome"),]
Simul4=Simul[(Simul$u==1e-08 & Simul$h==0.2 & Simul$s==-0.001 & Simul$InvSize==1000000 & Simul$N==10000 & Simul$Chromosome=="Ychromosome"),]



SimulSub=rbind(Simul1, Simul2, Simul3, Simul4) # Merge the data of the 3 focal simulations

SimulSub$Code=paste(SimulSub$N,SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$InvSize,SimulSub$Chromosome, SimulSub$Rep, sep="_") # Define a code identifying eavh simulation
summarySub=SimulSub %>% group_by(Code,N,u,InvSize, h, s, Rep, Chromosome) %>% summarise(maxFreq=last(Freq), maxGen=last(Gen), InitMutNumb=first(MeanMutInv))
summarySub=summarySub %>% group_by(N,u,InvSize, h, s, Chromosome) %>% mutate(GroupMaxGen=max(maxGen))
FixedSimul=summarySub[(summarySub$maxFreq>0.95 & summarySub$maxGen<6*summarySub$N) ,]$Code #Grep the code of the inversion that have fixed

LineToAdd2=SimulSub[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate ending states of fixed inversion.

for (i in FixedSimul){ #Do the same thing for fixed inversions
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  MaxGroupGen=summarySub[summarySub$Code==i,]$GroupMaxGen
  FalseEndGoodSimulSub=SimulSub[(SimulSub$Code==i & SimulSub$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimulSub$Gen=MaxGroupGen #Change it to 24991 (last recorded generation)
  FalseEndGoodSimulSub$Freq=1 #Set its frequency to 1
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub  #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
 b)
}

GoodSimulSubComplete=rbind(SimulSub,LineToAdd2) #Add these false simulation end to the simulation data.frame

GoodSimulSubComplete[GoodSimulSubComplete$Chromosome=="Autosome",]$Rep=100000+GoodSimulSubComplete[GoodSimulSubComplete$Chromosome=="Autosome",]$Rep
GoodSimulSubComplete$U=GoodSimulSubComplete$u*GoodSimulSubComplete$InvSize

GoodSimulSubComplete$Parameter <- with(GoodSimulSubComplete,
                                       paste0("bold(","N==", N,
                                              "~','~~s==", s,
                                              "~','~~h==", h,
                                              "~','~~U[L]==", U, ")")
)
options(scipen=0)
DataAllEnd=GoodSimulSubComplete %>% group_by(N,U, h, s, Rep, Chromosome) %>% summarise(MaxGen=max(Gen), MaxFreq=max(Freq), LastFreq=last(Freq))#Grep the end generation of each simulation
DataAllEnd$state="Lost"
DataAllEnd[(DataAllEnd$LastFreq==1),]$state="Fixed"
DataAllEnd[(DataAllEnd$MaxGen==6*DataAllEnd$N & DataAllEnd$MaxFreq<1),]$state="Segregating"
DataAllEnd = DataAllEnd %>% group_by(N,U, h, s, Chromosome) %>% mutate(GroupMaxGen=max(MaxGen))
SumEnd=DataAllEnd %>% group_by(N,U, h, s, Chromosome, GroupMaxGen) %>% count(state, sort = TRUE)  #Summarize the data for plotting summary number
SumEnd=SumEnd[SumEnd$state!="Lost",]
SumEnd$Pos=0.07
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chromosome=="Ychromosome"),]$Pos=1.0 #Define Position for plotting the text
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chromosome=="Autosome"),]$Pos=0.9

SumEnd$Parameter <- with(SumEnd,
                         paste0("bold(","N==", N,
                                "~','~~s==", s,
                                "~','~~h==", h,
                                "~','~~U[L]==", U, ")")
)

SumEnd$Parameter <- with(SumEnd,
                         paste0("bold(","N==", N,
                                "~','~~s==", s,
                                "~','~~h==", h,
                                "~','~~U[L]==", U,")")
)

Col2=c("#003049","#780000")

base=ggplot(GoodSimulSubComplete) 
PlotTrajDrift=base+rasterize(geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=1, size=0.2)+
  geom_text(data=SumEnd, aes(x=GroupMaxGen-500, y=Pos, label=paste0("n=",n), color=as.factor(Chromosome)), 
            vjust = -0.5, hjust = 0, size=4, show.legend = FALSE)+
  scale_color_manual("", values=rev(Col2), label=c("Autosomal","Y-linked"))+
  xlab("")+ylab("Inversion frequency")+
  ggtitle("Inversion trajectory")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face="bold"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    strip.text.x = element_text(size = 18, face="bold"),
    strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  #guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid( ~ Parameter, scale="free_x", labeller = labeller(Parameter=label_parsed))
PlotTrajDrift

#Mutation number


GoodSimulSubComplete=GoodSimulSubComplete %>% group_by(N,u,r,h,s,InvSize,Chromosome, U, Rep) %>% mutate(lastGen=max(Gen))
GoodSimulSubComplete=GoodSimulSubComplete %>% group_by(N,u,r,h,s,InvSize,Chromosome, U) %>% mutate(GroupMaxGen=max(lastGen))

df=GoodSimulSubComplete[!(GoodSimulSubComplete$Gen==GoodSimulSubComplete$lastGen),] #For plotting purpose, remove the last gen
# 1) Make a unique trajectory id (Chromosome + Rep)
df$traj <- interaction(df$Chromosome, df$Rep, drop = TRUE)
# 2) Randomize the order trajectories are drawn
df$traj <- factor(df$traj, levels = sample(levels(df$traj)))
# 3) Ensure points within each trajectory are ordered along x
df <- df[order(df$traj, df$Gen), ]
df$Log2MutNumber=log2(df$MeanMutInv/df$MeanMutNoInv)
df$Log2Fitness=log2(df$InvFit/df$NoInvFit)
df=df %>% group_by(N,u,r,h,s,InvSize,Chromosome, Rep, U,) %>% mutate(maxGen=max(Gen), MaxFreq=max(Freq), StartMutNumberNoInv=first(MeanMutNoInv),StartMutNumberInv=first(MeanMutInv))

### Figure 3C####
ColInv=c("#780000", "#003049")
ColNoInv=c("#c1121f", "#669bbc")

base=ggplot(df) 

PlotMutAccul=base+rasterize(geom_line(aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Mutation number 
  #geom_hline(yintercept = 0, linetype=1, size=0.2)+
  scale_color_manual("", values=ColInv, label=c("Y-linked"))+
  xlab("Generation")+ylab("Number of mutations")+
  #ThemeSobr+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    axis.line = element_line(colour = "grey"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16, face="bold"),
    strip.text.x =  element_blank(),
    legend.key.width = unit(1,"cm"),
    legend.key.height = unit(1, 'cm'),
    text = element_text(size=16),
    legend.position="none",
    axis.title.y = element_text(size=24, face="bold"),
    axis.title.x = element_blank()
  )+
  guides(color = guide_legend(override.aes = list( alpha=1)))+
  ggtitle("Mutation number in inverted segments")+
  facet_grid2( ~ Parameter, scales="free", independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotMutAccul


PlotMutAcculNoInv=base+rasterize(geom_line(aes(x=Gen, y=MeanMutNoInv, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Mutation number 
  #geom_hline(yintercept = 0, linetype=1, size=0.2)+
  scale_color_manual("", values=ColNoInv, label=c("Y-linked"))+
  xlab("Generation")+ylab("Number of mutations")+
  #ThemeSobr+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    axis.line = element_line(colour = "grey"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16, face="bold"),
    strip.text.x =  element_blank(),
    legend.key.width = unit(1,"cm"),
    legend.key.height = unit(1, 'cm'),
    text = element_text(size=16),
    legend.position="none",
    axis.title.y = element_text(size=24, face="bold"),
    axis.title.x = element_blank()
  )+
  ggtitle("Mutation number in non-inverted segments")+
  guides(color = guide_legend(override.aes = list( alpha=1)))+
  facet_grid2( ~ Parameter, scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotMutAcculNoInv



PlotTrajMutNumberLog=base+
  geom_hline(yintercept = 0, linetype="dashed")+
  rasterize(geom_line(aes(x=Gen, y=Log2MutNumber, group=traj), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=c('grey10'))+
  xlab("")+ylab("Log2 FC ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Inverted/Non-inverted mutation number ratio (Log2 Fold change)")+
  facet_grid2( ~ Parameter, scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))

base=ggplot(df[(df$Gen>0),]) 
PlotTrajFitnessInv=base+
  rasterize(geom_line(aes(x=Gen, y=InvFit, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColInv, label=c("Y-linked"))+
  xlab("")+ylab("fitness ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Fitness of inverted segments")+
  facet_grid2( ~ Parameter  , scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotTrajFitnessInv

PlotTrajFitnessNoInv=base+
  rasterize(geom_line(aes(x=Gen, y=NoInvFit, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColNoInv, label=c("Y-linked"))+
  xlab("")+ylab("fitness ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Fitness of non-inverted segments")+
  facet_grid2( ~ Parameter  , scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))


signed_sqrt_trans <- trans_new(
  "signed_sqrt",
  transform = function(x) sign(x) * sqrt(abs(x)),
  inverse   = function(x) sign(x) * (x^2)
)

PlotTrajFitnessLog=base+
  geom_hline(yintercept = 0, linetype="dashed")+
  rasterize(geom_line(aes(x=Gen, y=Log2Fitness, group=traj), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=c('grey10'))+
  xlab("")+ylab("log2 fold change ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  scale_y_continuous(trans = signed_sqrt_trans)+
  ggtitle("Inverted/Non-inverted fitness ratio (Log2 fold change)")+
  facet_grid2( ~ Parameter, scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))

PlotFreqNoInv=base+
  # geom_hline(yintercept = 0, linetype="dashed")+
  rasterize(geom_line(aes(x=Gen, y=FreqMutNoInv, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColNoInv, label=c("Y-linked"))+
  xlab("Generation")+ylab("Frequency ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_text(size = 24, face="bold"))+
  ggtitle("Mean frequency of mutations carried by the non-inverted segments ")+
  facet_grid2( ~ Parameter  , scale="free", independent = "y", labeller = labeller(Parameter=label_parsed))

PlotFreqInv=base+
  rasterize(geom_line(aes(x=Gen, y=FreqMutInv, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColInv, label=c("Y-linked"))+
  xlab("")+ylab("Frequency ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size = 24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Mean frequency of mutations carried by the inverted segments ")+
  facet_grid2( ~ Parameter  , scale="free",  independent = "y", labeller = labeller(Parameter=label_parsed))

PlotMerged=PlotTrajDrift / PlotMutAccul / PlotMutAcculNoInv / PlotTrajMutNumberLog / PlotTrajFitnessInv / PlotTrajFitnessNoInv / PlotTrajFitnessLog / PlotFreqInv / PlotFreqNoInv + plot_annotation(tag_levels = 'a')
#PlotMerged
save_plot("~/Project/MutationShelteringV2/Revision1/Plots/SubSet2_Y.png", PlotMerged, ncol=4, nrow=9)

###############################################@

# Figure S3 #
Simul1=Simul[(Simul$u==1e-08 & Simul$h==0.2 & Simul$s==-0.001 & Simul$InvSize==5000000 & Simul$N==1000 & Simul$Chromosome=="Autosome"),]
Simul2=Simul[(Simul$u==1e-08 & Simul$h==0.1 & Simul$s==-0.01 & Simul$InvSize==1000000 & Simul$N==1000 & Simul$Chromosome=="Autosome"),]
Simul3=Simul[(Simul$u==1e-08 & Simul$h==0.01 & Simul$s==-0.1 & Simul$InvSize==5000000 & Simul$N==10000 & Simul$Chromosome=="Autosome"),]
Simul4=Simul[(Simul$u==1e-08 & Simul$h==0.2 & Simul$s==-0.001 & Simul$InvSize==1000000 & Simul$N==10000 & Simul$Chromosome=="Autosome"),]

SimulSub=rbind(Simul1, Simul2, Simul3, Simul4) # Merge the data of the 3 focal simulations

SimulSub$Code=paste(SimulSub$N,SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$InvSize,SimulSub$Chromosome, SimulSub$Rep, sep="_") # Define a code identifying eavh simulation
summarySub=SimulSub %>% group_by(Code,N,u,InvSize, h, s, Rep, Chromosome) %>% summarise(maxFreq=last(Freq), maxGen=last(Gen), InitMutNumb=first(MeanMutInv))
summarySub=summarySub %>% group_by(N,u,InvSize, h, s, Chromosome) %>% mutate(GroupMaxGen=max(maxGen))
FixedSimul=summarySub[(summarySub$maxFreq>0.95 & summarySub$maxGen<6*summarySub$N) ,]$Code #Grep the code of the inversion that have fixed

LineToAdd2=SimulSub[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate ending states of fixed inversion.

for (i in FixedSimul){ #Do the same thing for fixed inversions
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  MaxGroupGen=summarySub[summarySub$Code==i,]$GroupMaxGen
  FalseEndGoodSimulSub=SimulSub[(SimulSub$Code==i & SimulSub$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimulSub$Gen=MaxGroupGen #Change it to 24991 (last recorded generation)
  FalseEndGoodSimulSub$Freq=1 #Set its frequency to 1
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub  #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
  #LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub)
}

GoodSimulSubComplete=rbind(SimulSub,LineToAdd2) #Add these false simulation end to the simulation data.frame
#GoodSimulSubComplete=SimulSub
GoodSimulSubComplete[GoodSimulSubComplete$Chromosome=="Autosome",]$Rep=100000+GoodSimulSubComplete[GoodSimulSubComplete$Chromosome=="Autosome",]$Rep
GoodSimulSubComplete$U=GoodSimulSubComplete$u*GoodSimulSubComplete$InvSize

GoodSimulSubComplete$Parameter <- with(GoodSimulSubComplete,
                                       paste0("bold(","N==", N,
                                              "~','~~s==", s,
                                              "~','~~h==", h,
                                              "~','~~U[L]==", U, ")")
)
options(scipen=0)

DataAllEnd=GoodSimulSubComplete %>% group_by(N,U, h, s, Rep, Chromosome) %>% summarise(MaxGen=max(Gen), MaxFreq=max(Freq), LastFreq=last(Freq))#Grep the end generation of each simulation
DataAllEnd$state="Lost"
DataAllEnd[(DataAllEnd$LastFreq==1),]$state="Fixed"
DataAllEnd[(DataAllEnd$MaxGen==6*DataAllEnd$N & DataAllEnd$MaxFreq<1),]$state="Segregating"
DataAllEnd = DataAllEnd %>% group_by(N,U, h, s, Chromosome) %>% mutate(GroupMaxGen=max(MaxGen))
SumEnd=DataAllEnd %>% group_by(N,U, h, s, Chromosome, GroupMaxGen) %>% count(state, sort = TRUE)  #Summarize the data for plotting summary number
SumEnd=SumEnd[SumEnd$state!="Lost",]
SumEnd$Pos=0.07
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chromosome=="Ychromosome"),]$Pos=1.0 #Define Position for plotting the text
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chromosome=="Autosome"),]$Pos=0.9




SumEnd$Parameter <- with(SumEnd,
                         paste0("bold(","N==", N,
                                "~','~~s==", s,
                                "~','~~h==", h,
                                "~','~~U[L]==", U,")")
)


Col2=c("#780000", "#003049")
base=ggplot(GoodSimulSubComplete) 
PlotTrajDrift=base+rasterize(geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=1, size=0.2)+
  geom_text(data=SumEnd, aes(x=GroupMaxGen-500, y=Pos, label=paste0("n=",n), color=as.factor(Chromosome)), 
            vjust = -0.5, hjust = 0, size=4, show.legend = FALSE)+
  scale_color_manual("", values=rev(Col2), label=c("Autosomal","Y-linked"))+
  xlab("")+ylab("Inversion frequency")+
  ggtitle("Inversion trajectory")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face="bold"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    strip.text.x = element_text(size = 18, face="bold"),
    strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  #guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid( ~ Parameter, scale="free_x", labeller = labeller(Parameter=label_parsed))

GoodSimulSubComplete=GoodSimulSubComplete %>% group_by(N,u,r,h,s,InvSize,Chromosome, U, Rep) %>% mutate(lastGen=max(Gen))
GoodSimulSubComplete=GoodSimulSubComplete %>% group_by(N,u,r,h,s,InvSize,Chromosome, U) %>% mutate(GroupMaxGen=max(lastGen))

df=GoodSimulSubComplete[!(GoodSimulSubComplete$Gen==GoodSimulSubComplete$lastGen),] #For plotting purpose, remove the last gen
#df <- GoodSimulSubComplete[!(GoodSimulSubComplete$Gen==6*GoodSimulSubComplete$N & GoodSimulSubComplete$Freq==1),] # replace with your plotting df
# 1) Make a unique trajectory id (Chromosome + Rep)
df$traj <- interaction(df$Chromosome, df$Rep, drop = TRUE)
# 2) Randomize the order trajectories are drawn
df$traj <- factor(df$traj, levels = sample(levels(df$traj)))
# 3) Ensure points within each trajectory are ordered along x
df <- df[order(df$traj, df$Gen), ]
df$Log2MutNumber=log2(df$MeanMutInv/df$MeanMutNoInv)
df$Log2Fitness=log2(df$InvFit/df$NoInvFit)
df=df %>% group_by(N,u,r,h,s,InvSize,Chromosome, Rep, U,) %>% mutate(maxGen=max(Gen), MaxFreq=max(Freq), StartMutNumberNoInv=first(MeanMutNoInv),StartMutNumberInv=first(MeanMutInv))

ColInv=c( "#003049")
ColNoInv=c("#669bbc")
base=ggplot(df) 

PlotMutAccul=base+rasterize(geom_line(aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Mutation number 
  scale_color_manual("", values=ColInv, label=c("Y-linked"))+
  xlab("Generation")+ylab("Number of mutations")+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    axis.line = element_line(colour = "grey"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16, face="bold"),
    strip.text.x =  element_blank(),
    legend.key.width = unit(1,"cm"),
    legend.key.height = unit(1, 'cm'),
    text = element_text(size=16),
    legend.position="none",
    axis.title.y = element_text(size=24, face="bold"),
    axis.title.x = element_blank()
  )+
  guides(color = guide_legend(override.aes = list( alpha=1)))+
  ggtitle("Mutation number in inverted segments")+
  facet_grid2( ~ Parameter, scales="free", independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotMutAccul


PlotMutAcculNoInv=base+rasterize(geom_line(aes(x=Gen, y=MeanMutNoInv, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Mutation number 
  scale_color_manual("", values=ColNoInv, label=c("Y-linked"))+
  xlab("Generation")+ylab("Number of mutations")+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    axis.line = element_line(colour = "grey"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16, face="bold"),
    strip.text.x =  element_blank(),
    legend.key.width = unit(1,"cm"),
    legend.key.height = unit(1, 'cm'),
    text = element_text(size=16),
    legend.position="none",
    axis.title.y = element_text(size=24, face="bold"),
    axis.title.x = element_blank()
  )+
  ggtitle("Mutation number in non-inverted segments")+
  guides(color = guide_legend(override.aes = list( alpha=1)))+
  facet_grid2( ~ Parameter, scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotMutAcculNoInv



PlotTrajMutNumberLog=base+
  geom_hline(yintercept = 0, linetype="dashed")+
  rasterize(geom_line(aes(x=Gen, y=Log2MutNumber, group=traj), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=c('grey10'))+
  xlab("")+ylab("Log2 FC ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Inverted/Non-inverted mutation number ratio (Log2 Fold change)")+
  facet_grid2( ~ Parameter, scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))


base=ggplot(df[(df$Gen>0),]) 
PlotTrajFitnessInv=base+
  rasterize(geom_line(aes(x=Gen, y=InvFit, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColInv, label=c("Y-linked"))+
  xlab("")+ylab("fitness ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    #strip.text = element_text(face="bold"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    #strip.text.x = element_text(size = 14, face="bold"),
    #strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Fitness of inverted segments")+
  facet_grid2( ~ Parameter  , scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotTrajFitnessInv

PlotTrajFitnessNoInv=base+
  rasterize(geom_line(aes(x=Gen, y=NoInvFit, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColNoInv, label=c("Y-linked"))+
  xlab("")+ylab("fitness ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Fitness of non-inverted segments")+
  facet_grid2( ~ Parameter  , scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))


signed_sqrt_trans <- trans_new(
  "signed_sqrt",
  transform = function(x) sign(x) * sqrt(abs(x)),
  inverse   = function(x) sign(x) * (x^2)
)


PlotTrajFitnessLog=base+
  geom_hline(yintercept = 0, linetype="dashed")+
  rasterize(geom_line(aes(x=Gen, y=Log2Fitness, group=traj), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=c('grey10'))+
  xlab("")+ylab("log2 fold change ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_blank())+
  scale_y_continuous(trans = signed_sqrt_trans)+
  ggtitle("Inverted/Non-inverted fitness ratio (Log2 fold change)")+
  facet_grid2( ~ Parameter, scale="free",    independent = "y", labeller = labeller(Parameter=label_parsed))


###

PlotFreqNoInv=base+
  rasterize(geom_line(aes(x=Gen, y=FreqMutNoInv, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColNoInv, label=c("Y-linked"))+
  xlab("Generation")+ylab("Frequency ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size=24, face="bold"),
    axis.title.x = element_text(size = 24, face="bold"))+
  ggtitle("Mean frequency of mutations carried by the non-inverted segments ")+
  facet_grid2( ~ Parameter  , scale="free", independent = "y", labeller = labeller(Parameter=label_parsed))
#PlotFreqNoInv

PlotFreqInv=base+
  rasterize(geom_line(aes(x=Gen, y=FreqMutInv, group=traj, color=as.factor(Chromosome)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  scale_color_manual("", values=ColInv, label=c("Y-linked"))+
  xlab("")+ylab("Frequency ")+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=24, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size = 24, face="bold"),
    axis.title.x = element_blank())+
  ggtitle("Mean frequency of mutations carried by the inverted segments ")+
  facet_grid2( ~ Parameter  , scale="free",  independent = "y", labeller = labeller(Parameter=label_parsed))

PlotMerged=PlotTrajDrift / PlotMutAccul / PlotMutAcculNoInv / PlotTrajMutNumberLog / PlotTrajFitnessInv / PlotTrajFitnessNoInv / PlotTrajFitnessLog / PlotFreqInv / PlotFreqNoInv + plot_annotation(tag_levels = 'a')
save_plot("~/Project/MutationShelteringV2/Revision1/Plots/SubSet2_Autosome.png", PlotMerged, ncol=4, nrow=9)


#########################

### Figure S4 ###
# See file Figure3_V070526.R Figure3_and_S4_V070526.R 

#########################

### Figure S5 ###

Simul=read.table(paste("~/Project/MutationShelteringV2/Revision1/Output/IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat.24Col.txt",sep=""), stringsAsFactors = F) # Data for N=10000
colnames(Simul)=c("N", "u", "r", "h", "s", "sInv", "Gen", "StartInv", "EndInv", "Rep" ,"Chromosome",
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")

Simul$Gen=Simul$Gen - 6*Simul$N
Simul[Simul$Chromosome=="SexChrom",]$Freq=Simul[Simul$Chromosome=="SexChrom",]$Freq * 4 # Define Y inversion frequency
Simul$InvSize=Simul$EndInv - Simul$StartInv #Inversion size
Simul[Simul$Chromosome=="SexChrom",]$Chromosome="Ychromosome"

Summary2=Simul %>% group_by(N, u, r, h, s, InvSize, Chromosome, Rep) %>% summarise(maxFreq=last(Freq), maxGen=last(Gen), InitMutNumb=first(MeanMutInv), InitMutNumbNoInv=first(MeanMutNoInv)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)
#Summary2=summarySub2
Summary2$State="LostEarly"
Summary2[Summary2$maxGen==6*Summary2$N,]$State="Segregating"
Summary2[Summary2$maxFreq>0.95,]$State="Fixed"
#Summary2$log2FC=log2(Summary2$InitMutNumb/Summary2$InitMutNumbNoInv)
Summary2$sNorm=Summary2$s*Summary2$N
Summary2$FixCode=0
Summary2[Summary2$State=="Fixed",]$FixCode=1

ParamSummary2=Summary2 %>% group_by(N, u, r, h, s, sNorm, InvSize, Chromosome) %>% summarise(ProbFix=mean(FixCode), InitMut=mean(InitMutNumbNoInv))
Col7=scales::viridis_pal(begin=0.1, end=0.9, option="C")(7)
options(scipen = 999)
ParamSummary2$U=ParamSummary2$u * ParamSummary2$InvSize
PlotMutSegFix=ggplot(ParamSummary2[(ParamSummary2$sNorm %in% c(-1,-10,-100) & ParamSummary2$U %in% c(0.001, 0.01, 0.05) & ParamSummary2$Chromosome=="Ychromosome"),], aes(x=InitMut, y=ProbFix, color=as.factor(h), shape=as.factor(sNorm), size=as.factor(U)))+
  geom_point(alpha=0.8)+
  scale_x_log10("Average number of derived mutations per invividuals \n in the inversion region (log scale)")+
  scale_y_log10("Fixation probability (log10 scale)")+
  scale_color_manual("h",values=Col7)+
  scale_size_manual(expression(U[L]), values=c(2,3,4))+
  scale_shape_manual("Ns", values=c(15,16,17))+ThemeSobr+
  theme(panel.grid.major=element_line(color="grey95"))
PlotMutSegFix
save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/FixProbaMutNumber.pdf"), PlotMutSegFix, ncol=2, nrow=2) 

########################################

### Figure S6 ###

Simul=read.table(paste("~/Project/MutationShelteringV2/Revision1/Output/IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat.24Col.txt",sep=""), stringsAsFactors = F) # Data for N=10000
colnames(Simul)=c("N", "u", "r", "h", "s", "sInv", "Gen", "StartInv", "EndInv", "Rep" ,"Chromosome",
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")

Simul$Gen=Simul$Gen - 6*Simul$N
Simul[Simul$Chromosome=="SexChrom",]$Freq=Simul[Simul$Chromosome=="SexChrom",]$Freq * 4 # Define Y inversion frequency
Simul$InvSize=Simul$EndInv - Simul$StartInv #Inversion size
Simul[Simul$Chromosome=="SexChrom",]$Chromosome="Ychromosome"

Summary2=Simul %>% group_by(N, u, r, h, s, InvSize, Chromosome, Rep) %>% summarise(maxFreq=last(Freq), maxGen=last(Gen), InitMutNumb=first(MeanMutInv), InitMutNumbNoInv=first(MeanMutNoInv)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)
#Summary2=summarySub2
Summary2$State="LostEarly"
Summary2[Summary2$maxGen==6*Summary2$N,]$State="Segregating"
Summary2[Summary2$maxFreq>0.95,]$State="Fixed"
#Summary2$log2FC=log2(Summary2$InitMutNumb/Summary2$InitMutNumbNoInv)
Summary2$sNorm=Summary2$s*Summary2$N

write.table(Summary2, "~/Project/MutationShelteringV2/Revision1/Output/SummaryEnd.tsv",sep="\t", quote=F, row.names=F ) # Save dataset, to save time, just in case

DataAllEnd=read.table("~/Project/MutationShelteringV2/Revision1/Output/SummaryEnd.tsv",sep="\t",header=T )

DataAllEndFixed=DataAllEnd[(DataAllEnd$State=="Fixed"),] 
DataAllEndFixed$U=DataAllEndFixed$InvSize * DataAllEndFixed$u
DataAllEndFixed$Parameter <- with(DataAllEndFixed,
                                   paste0("bold(","N==", N,
                                          "~','~~U[L]==", U, ")")
)
Col3=scales::viridis_pal(begin=0.2, end=0.8, option="C")(3)
base=ggplot(DataAllEndFixed[(DataAllEndFixed$U %in% c(0.01, 0.001, 0.05) & DataAllEndFixed$h>0 & DataAllEndFixed$sNorm %in% c(-1, -10, -100)),])
PlotFixTime=base+
  #geom_hline(aes(yintercept=1), linetype="dashed", color="grey")+
  geom_boxplot(aes(fill=as.factor(sNorm), y=maxGen, x=as.factor(h) ), outliers = F)+
  facet_grid(
    Chromosome ~ Parameter,
    labeller = label_parsed
  )+
  scale_fill_manual(expression(paste("Selection coefficient (",italic("Ns"), ")")), values=Col3)+
  ThemeSobr+
  theme(panel.border = element_blank(),
        # legend.position = c(0.5,1.1),
        # legend.direction = "horizontal",
        panel.grid.major = element_line(colour="grey97"),
        legend.text = element_text(face="bold", size=18),
        legend.title = element_text(face="bold", size=18),
        strip.text = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=18),
        axis.text = element_text(size=16),
        legend.background = element_blank(),
        plot.title=element_text(size=11, face="bold",hjust=0.5),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        plot.margin = margin(29, 3, 3, 5, "pt"))+
  ylab("Fixation time (Generation, log10 scale)")+
  scale_y_log10()+
  xlab(expression(paste(bold("Dominance coefficient ("),bolditalic("h"), bold(")"))))
#labs(title="Relative number of mutations in fixed inversions")
PlotFixTime
save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/PlotFixTime.pdf"), PlotFixTime, nrow=2, ncol=6, base_aspect_ratio = 1)    


###########################################################

### Figure S7 ###
#See file Figure4_and_S7_V120326.R

######################

### Figure S8 ; Recombination Suppressor with PHL ###
Col=scales::viridis_pal(begin=0, end=0.6, option="A")(2)
Simul=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/RecombinationSuppressorEvolution_n=2Mb_FigS15.txt",sep=""), stringsAsFactors = F) #Evolution of recombination suppressors (suppress recombination in heterozygote and in homozygote, unlike inversions)
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
#Simul=subset(Simul, Simul$Gen %% 10 == 0) #Keep only every 10 generation for computing purpose
Simul$InvSize=Simul$FinInv - Simul$DebInv #Size of the region affected (n)
SimulSub=Simul #Copy, just in case
SimulSub$Link="Linked" #Define linkage
SimulSub[SimulSub$DebInv > 10000000,]$Link="Unlinked" 
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv, SimulSub$Rep, sep="_") #Define Code
FocS=c(-0.01,-0.05,-0.1) #Value of s to focus on 
SimulSub=SimulSub[(SimulSub$s %in% FocS),] 
#SimulSub20G=unique(SimulSub[SimulSub$Gen>15020,]$Code) #Keep only inversion not lost during the first 20 generation
#SimulSub=SimulSub[SimulSub$Code %in% SimulSub20G,]
summarySub=SimulSub %>% group_by(Code,h, s, InvSize, u, N,Link) %>% summarise(LastFreq=last(Freq), maxGen=last(Gen)) # As before, recreate end state
SummaryAll=summarySub %>% group_by(h, s, InvSize, u, N, Link) %>% summarise(n=n()) # As before, recreate end state
LostSimulSub=summarySub[summarySub$maxGen<25000,]$Code # Inversion that have been lost or fixed
GoodSimulSub=SimulSub #Keep only non-lost Inversion
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15002),]#For lost inversion, grep their initial state
FalseEndGoodSimulSub$Gen=25000 # Modify their initial state
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
#FixedSimul=summarySub[summarySub$LastFreq>0.95,]$Code #Grep the code of the inversion that have fixed
#FalseEndGoodSimulSub[FalseEndGoodSimulSub$Code %in% FixedSimul,]$Freq=1.0 #For all inversion that have fixed, define their end frequency as 1.0
GoodSimulSub=rbind(GoodSimulSub,FalseEndGoodSimulSub) #Concatenate the dataset
DataAll=GoodSimulSub
DataAll$s=-DataAll$s
DataAllEnd=subset(DataAll, DataAll$Gen==25000)
DataAllEnd$state="Lost"
DataAllEnd[DataAllEnd$Freq==0.5,]$state="Segregate"
#SummaryAllEnd=GoodSimulSub %>% group_by(h, s, InvSize, u, N, Link) %>% summarise(n=n())
SumEnd=DataAllEnd %>% count(h, s, InvSize, u, N, Link, state, sort = TRUE) 
SumEnd$Pos=0.0
SumEnd[SumEnd$state=="Segregate",]$Pos=0.5

DataAll$Gen=DataAll$Gen-15000 #Suppress the burn-in generation.
themeInvFreq=theme(legend.position = c(0.50,0.98), #Change a bit the theme.
                   legend.direction = "horizontal",
                   panel.border = element_blank(),  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   text = element_text(size=25),
                   axis.line = element_line(colour = "grey"),
                   strip.text.x  = element_blank(),
                   strip.text.y  = element_text(size=25),
                   strip.background.x = element_blank() )

base=ggplot(DataAll[(DataAll$h==0.1 & DataAll$InvSize==2000000),])
Plot0.1_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Recombination suppressor frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.1 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.1")
#Plot0.1_2Mb

base=ggplot(DataAll[(DataAll$h==0.01 & DataAll$InvSize==2000000),])
Plot0.01_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Recombination suppressor frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.01 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.01")

base=ggplot(DataAll[(DataAll$h==0.001 & DataAll$InvSize==2000000),])
Plot0.001_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Recombination suppressor requency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.001 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.001")

Plot2Mb=Plot0.001_2Mb / Plot0.01_2Mb / Plot0.1_2Mb + plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(size = 16, face="bold"))

save_plot("~/Project/MutationShelteringV2/Revision1/Plots/FigS7_RecombSuppressor.png", Plot2Mb, ncol=2, nrow=9)
save_plot("~/Project/MutationShelteringV2/Revision1/Plots/FigS7_RecombSuppressor.pdf", Plot2Mb, ncol=4, nrow=6)

##################################

### Figure S9 ###
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), s2=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(), FYN=double(),FYI=double(), FY=double(),
                      W=double(),D=double(),q=double())

write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_Over.txt", append = F, quote=F, row.names = F)
P=0.8
R=0.0
for (s2 in c(0.0, 0.5, 0.9, 0.95, 0.99,1))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.05, 0.1, 0.5))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), s2=double(), P=double(),
                            FXN=double(),FXI=double(),FX=double(), FYN=double(),FYI=double(), FY=double(),
                            W=double(),D=double(),q=double())
      
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXN=0.49
      FXI=0.01
      FYN=0.50
      FYI=0.0
      FY=FYN+FYI
      FX=FXN+FXI
      W=WII*((1-s2)*(FXI*FXI + FYI*FYI) + 2*FXI*FYI) + 
        WNN*((1-s2)*(FXN*FXN + FYN*FYN) + 2*FXN*FYN) +
        2*WNI*((1-s2)*(FXI*FXN + FYI*FYN) + FXI*FYN + FXN*FYI)
      D=FXI*FYN - FXN*FYI
      FirstGen=c(1,h,s,u,n,s2,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        FXI=(TableLinkR$FXI[time-1]*(WII*((1-s2)*TableLinkR$FXI[time-1] + TableLinkR$FYI[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FXN[time-1] + TableLinkR$FYN[time-1])) - 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FYI=(TableLinkR$FYI[time-1]*(WII*((1-s2)*TableLinkR$FYI[time-1] + TableLinkR$FXI[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FYN[time-1] + TableLinkR$FXN[time-1])) + 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FXN=(TableLinkR$FXN[time-1]*(WNN*((1-s2)*TableLinkR$FXN[time-1] + TableLinkR$FYN[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FXI[time-1] + TableLinkR$FYI[time-1])) + 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FYN=(TableLinkR$FYN[time-1]*(WNN*((1-s2)*TableLinkR$FYN[time-1] + TableLinkR$FXN[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FYI[time-1] + TableLinkR$FXI[time-1])) - 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FY=FYN+FYI
        FX=FXN+FXI
        W=WII*((1-s2)*(FXI*FXI + FYI*FYI) + 2*FXI*FYI) + 
          WNN*((1-s2)*(FXN*FXN + FYN*FYN) + 2*FXN*FYN) +
          2*WNI*((1-s2)*(FXI*FXN + FYI*FYN) + FXI*FYN + FXN*FYI)
        D=FXI*FYN - FXN*FYI
        Gen=c(time,h,s,u,n,s2,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_Over.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8)(6)
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_Over.txt", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1)
Data=Data[Data$s %in% FocS,]
base=ggplot(Data)
Plot0.8=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=FXI, color=fct_rev(fct_infreq(as.factor(s2)))), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("s2=", values=Col)+
  theme(panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=25),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5), limits=c(0,0.6))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))#+

save_plot("FigS9.png", Plot0.8, ncol = 3, nrow=4)
save_plot("FigS9.pdf", Plot0.8, ncol = 3, nrow=4)


### Figure S10 ###
#### With distribution of s values ###

SimulLamb=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/MutationLambdaDistribution_InversionTrajectories_FigS16.txt",sep=""), stringsAsFactors = F) #File available from https://doi.org/10.6084/m9.figshare.19961033
colnames(SimulLamb)=c("N", "u", "r", "s", "Gen", "DebInv", "FinInv", "Rep", 
                      "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                      "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                      "InvFit", "NoInvFit","Freq","Chromosome")

SimulLamb$Position="Y"
SimulLamb[SimulLamb$DebInv>10000000,]$Position="Autosome"
SimulLamb[SimulLamb$Position=="Y",]$Freq=SimulLamb[SimulLamb$Position=="Y",]$Freq * 4
SimulLamb$InvSize=SimulLamb$FinInv - SimulLamb$DebInv
SimulLambSub=SimulLamb
summarySub_BLamb=SimulLambSub %>% group_by(N,u,r,s,InvSize,Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv), MinSegMut=min(MeanMutNoInv))
summarySub_BLamb$State="LostEarly"
summarySub_BLamb[summarySub_BLamb$maxGen>15020,]$State="LostLate"
summarySub_BLamb[summarySub_BLamb$maxGen==24991,]$State="Segregating"
summarySub_BLamb[summarySub_BLamb$maxFreq>0.95,]$State="Fixed"
summarySub_BLamb$StateCode=1
summarySub_BLamb[summarySub_BLamb$State=="LostEarly",]$StateCode=0
summarySub_BLamb[summarySub_BLamb$State=="LostLate",]$StateCode=0
summarySub_BLamb[summarySub_BLamb$State=="Segregating",]$StateCode=0
DataSummaryLamb=summarySub_BLamb %>% group_by(N,u,r,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut))

options(scipen=0)
Col7=scales::viridis_pal(begin=0, end=0.9, option="A")(7)
options(scipen=999)
DataSubLamb=DataSummaryLamb[(DataSummaryLamb$s<0),]
DataSubLamb$s=-DataSubLamb$s
DataSubLamb$U=DataSubLamb$u*DataSubLamb$InvSize
DataSubLamb=DataSubLamb[!(DataSubLamb$U==0.005 & DataSubLamb$u==1e-8),] #Two combination of u and InvSize create a region with U=0.005. For plotting purpose, just remove one on these combinations.
DataSubLamb$NeutralFixProb=1/(2*DataSubLamb$N)
DataSubLamb[DataSubLamb$Position=="Y",]$NeutralFixProb=1/(DataSubLamb[DataSubLamb$Position=="Y",]$N/2)
DataSubLamb$RelatProbSpread=DataSubLamb$ProbSpread/DataSubLamb$NeutralFixProb

Base=ggplot(DataSubLamb[DataSubLamb$s<0.5,], aes( y=RelatProbSpread))
PlotPropInvSpread_AllInvLamb=Base+geom_point(aes(x=s, color=as.factor(U)), size=5)+
  geom_line(aes(x=s, color=as.factor(U)), size=2)+
  scale_color_manual("Region-wide mutation rate (U)", values=Col7)+
  scale_x_log10(breaks=c(0.1, 0.01, 0.001))+
  ylab("Normalised inversion fixation probability")+
  xlab("mean(s)")+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    legend.spacing.y= unit(0, 'cm'),
    panel.spacing = unit(0.8, "lines"),
    text = element_text(size=18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18))+
  facet_grid( ~ Position)
PlotPropInvSpread_AllInvLamb
save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/FigS?_sDistrib.png"),PlotPropInvSpread_AllInvLamb, ncol=2.3, nrow=3) #Fig. S16
save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/FigS?_sDistrib.pdf"),PlotPropInvSpread_AllInvLamb, ncol=2.3, nrow=3) #Fig. S16
######################

### Figure S11###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/HaploDiplo_InvTrajectories_FigS22.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000

colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "FreqHaplo", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")
Simul$InvSize=Simul$FinInv - Simul$DebInv # Inversion size
Simul$Position="Sex-linked inversions"
Simul[Simul$DebInv==500000,]$Position="Autosomal inversions"
Simul$Code=paste(Simul$N,Simul$u,Simul$r,Simul$h,Simul$s,Simul$InvSize,Simul$FreqHaplo,Simul$Position, Simul$Rep, sep="_")
Summary=Simul %>% group_by(N,u,r,h,s,InvSize,FreqHaplo, Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen)) # For each simulation,  grep its last generation (when the inversion was lost or fixed) and its maximum frequency
Summary$State="LostEarly" #Define an state defining inversion
Summary[Summary$maxGen>15020,]$State="LostLate" #Inversion lost after 20 generation or more
Summary[Summary$maxGen==24991,]$State="Segregating" #Inversion still segregating at simulation end
Summary[Summary$maxFreq==0.5,]$State="Fixed" #Inversion that reached above 0.95 frequency are considered fixed (for computation purpose, simulation stop when inversion fix,so sometime we do not observe inversion at 1.0)
Summary$StateCode=1 #For estimating the proportion of inversion fixed, note as 1 inversion fixed and 0 otherwise
Summary[Summary$State=="LostEarly",]$StateCode=0
Summary[Summary$State=="LostLate",]$StateCode=0
Summary[Summary$State=="Segregating",]$StateCode=0

SumNoLostEarly=subset(Summary, Summary$State!="LostEarly") #Remove inversion that were lost in fewer than 20 generations
DataSummary=SumNoLostEarly %>% group_by(N,u,r,h,s,InvSize, FreqHaplo, Position) %>% summarise(ProbSpread=mean(StateCode)) # For each set of parameter, compute the fraction of mutation fixed (only for not mutation-free inversion)
options(scipen=0) #Scientific notation
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(4) #Define color

Base=ggplot(DataSummary, aes( y=ProbSpread)) #Plot the result
PlotPropInvSpread=Base+geom_point(aes(x=as.factor(h), color=as.factor(FreqHaplo)), size=3, alpha=0.7)+
  scale_color_manual("Life cycle", values=Col, labels=c("Haplodiplontic (1/2 haploid)","Haplodiplontic (1/3 haploid)","Haplodiplontic (1/10 haploid)", "Fully diploid"))+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab("Fraction of inversions fixed \n after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  facet_grid(.~Position)+
  theme(panel.border = element_blank(), 
        legend.text = element_text(face="bold"),
        legend.key = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        legend.background = element_blank(),
        axis.title = element_text(face="bold"),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=4))
PlotPropInvSpread

save_plot("FigS11.png", PlotPropInvSpread, ncol=2, base_aspect_ratio = 1)
save_plot("FigS11.pdf", PlotPropInvSpread, ncol=2, base_aspect_ratio = 1)

####################


### Figure S12 ###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/N=1000_2MbChromosomeFusion_Trajectory_FigS18.txt",sep=""), stringsAsFactors = F) #File available from https://doi.org/10.6084/m9.figshare.19961033
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")

Simul[Simul$Chromosome=="Y",]$Freq=Simul[Simul$Chromosome=="Y",]$Freq * 4 #frequency of Y inversions in the population of Y chromosome, not the overall frequency
Simul[Simul$Chromosome=="X",]$Freq=Simul[Simul$Chromosome=="X",]$Freq * (1/0.75) #frequency of Y inversions in the population of Y chromosome, not the overall frequency

Simul$InvSize=Simul$FinInv - Simul$DebInv # Inversion size
Simul$Gen = Simul$Gen - 15000
Col=scales::viridis_pal(begin=0.2, end=0.8, option="A")(2)
base=ggplot(Simul) #Plot the data
PlotA=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chromosome)), size=0.5, alpha=0.3)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, label=c("X-autosome fusion","Y-autosome fusion"))+
  xlab("Generation")+ylab("Fused chromosome frequency")+
  ThemeSobr+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold", size=16),
        legend.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title = element_text(face="bold", size=16),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
        strip.placement = "outside",
        legend.key=element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
  )+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("h=", h) ~ Chromosome)

save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS18.pdf"),PlotA, nrow=3, ncol=2) #Fig. S18
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS18.png"),PlotA, nrow=3, ncol=2) #Fig. S18


### Figure S13 ###
# See file Figure5_and_S13.R

###############@

### Figure S14:  Burn-in stat ##
data=read.table("~/Project/MutationShelteringV2/Revision1/InitialState/slim_g6N_XYsyst_5M_Autosome_N1000-10000_r1.0e-07_uAll_NoDistrib_sAll_hAll_Rep1.MutFreq.txt", stringsAsFactors = F) #Stat computed during the burn in of each simulations
colnames(data)=c("N", "mu","r","s","generation", "Chromosome", "Rep", "meanNbMut", "Nmut1", "MeanFreq1", "h")
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(7)#Color Palette 
options(scipen=0) #Non-scientific notation


ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12),
  axis.line = element_line(colour = "grey"),
  legend.spacing.y= unit(0, 'cm'),
  panel.spacing = unit(0.8, "lines"),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 9)
)
base1000=ggplot(data[data$N==1000,]) #Plot number of mutation
nMut1000=base1000+geom_line(aes(x=generation, y=Nmut1, color=as.factor(h)))+ 
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Number of mutations genome-wide")+
  xlab("Generation (burn-in)")+
  ggtitle("N=1000")+
  ThemeSobr

FreqMut1000=base1000+geom_line(aes(x=generation, y=MeanFreq1, color=as.factor(h)))+
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Mean frequency of mutations")+
  xlab("Generation (burn-in)")+
  ggtitle("N=1000")+
  ThemeSobr

base10000=ggplot(data[data$N==10000,]) #Plot number of mutation
nMut10000=base1000+geom_line(aes(x=generation, y=Nmut1, color=as.factor(h)))+ 
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Number of mutations genome-wide")+
  xlab("Generation (burn-in)")+
  ggtitle("N=10000")+
  ThemeSobr

FreqMut10000=base10000+geom_line(aes(x=generation, y=MeanFreq1, color=as.factor(h)))+
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Mean frequency of mutations")+
  xlab("Generation (burn-in)")+
  ggtitle("N=10000")+
  ThemeSobr

PlotMerg= nMut1000 / FreqMut1000 / nMut10000 / FreqMut10000 + plot_annotation(tag_levels = list(c("a", "b", "c", "d"))) &
  theme(plot.tag = element_text(size = 16, face="bold"))
  
PlotMerg
save_plot("~/Project/MutationShelteringV2/Revision1/Plots/FigS?_BurnInStat.png", nrow = 8, ncol=3,PlotMerg)
save_plot("~/Project/MutationShelteringV2/Revision1/Plots/FigS?_BurnInStat.pdf", nrow = 8, ncol=3,PlotMerg)


