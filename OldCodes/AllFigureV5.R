#### Papier Sheltering, V3 (2025) 

library(cowplot)
library(ggplot2)
library(tidyverse)
library(viridis)
library(directlabels)
options(scipen = 0)
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


## Figure 2 #
# Panel A & B#
ResultHeatNonNeutral=data.frame(h=double(), s=double(), u=double(), U=double(), n=double(), q=double(), ### Table containing the result to display. 
                                nq=double(), TotProb=double(), CumulProbFix=double(), CumulSelectAdv=double(), FracFixedY=double(), FracFixedAuto=double())

u=2.2/3100000000 #Deleterious mutation rate in human.
for (U in c(0.001, 0.0001, 0.01, 0.05))
{
  n=as.integer(U/u)
  for (h in seq(0.001,0.5,0.0025))
  {
    for (s in  seq(0.001,0.1,0.001))
    {
      print(paste(c(U, h, s)))
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2)))) #Mutation frequencies
      nq=n*q #Average number of mutation of segment of size n.
      #B1=(q*n*h)/(1-h) #Threshold under which inversions fix on the autosome. (WNI > WNN)
      TotProb=0
      CumulProbNonNeutral=0
      CumulProbFix=0
      FracFixedY=0
      CumulSelectAdv=0
      FracFixedY=0 #Fraction of fixed inversion
      FracFixedAuto=0 #Fraction of fixed autosomal inversion
      for (m in seq(0,nq)) #here, consider less-loaded inversion. Replace by "m in seq(0,n)" to consider all inversions (Figure S9)
      {
        Pm=dbinom(m,n,q) #Probabilities of occuring inversion with m mutations
        #print(paste(c(m, n, q)))
        #print(Pm)
        if (m<nq) #If the inversion is less loaded... (which is always true here because "m in seq(0,nq-1)")
        {
          FracFixedY=FracFixedY+Pm
          WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
          WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
          WII=(1-s)^m #Fitness of individual homozygous for the inversion.
          DiffFitness=WNI-WNN
          if (WNI<WII)
          {
            FracFixedAuto=FracFixedAuto+Pm
          }
          AverageSelectAdv=Pm*DiffFitness
          ProbFix=Pm*2*DiffFitness
          CumulProbFix=CumulProbFix+ProbFix
          CumulSelectAdv=CumulSelectAdv+AverageSelectAdv
          TotProb=TotProb + Pm #Probabilities of less loaded inversion
        }
      }
      ResultHeatNonNeutral[nrow(ResultHeatNonNeutral)+1,]=c(h,s,u, U, n, q, n*q,TotProb, CumulProbFix, CumulSelectAdv, FracFixedY, FracFixedAuto) #Fill the table
    }
  }
}

write.table(ResultHeatNonNeutral,"~/Paper/ModelSexChrom_V2/Plots/DeterministicSimul.tsv", sep="\t", row.names = F, quote = F)
ResultHeatNonNeutral=read.table("~/Paper/ModelSexChrom_V2/Plots/DeterministicSimul.tsv", header=T, stringsAsFactors = F)

Subset=ResultHeatNonNeutral[( ResultHeatNonNeutral$h==0.1010 & ResultHeatNonNeutral$s==0.001 & ResultHeatNonNeutral$U==0.001 ),]
mean(ResxultHeatNonNeutral[(ResultHeatNonNeutral$U==0.005 ),]$CumulSelectAdv)
mean(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.001 ),]$CumulSelectAdv)

ResultHeatNonNeutral$s=-ResultHeatNonNeutral$s
options(scipen = 9)
base=ggplot(ResultHeatNonNeutral[ResultHeatNonNeutral$U %in% c(0.05, 0.01, 0.001),])
# Panel B
SelectAdv=base+geom_tile(aes(x=h, y=s, fill=CumulSelectAdv))+
  ylab(expression(paste("Selection coefficient (",italic("s"), ")")))+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))+
  # coord_cartesian(ylim=c(0,0.1), expand = F)+
  facet_grid(. ~ paste0("U=", U))+
  scale_fill_viridis(option = "B", "Average selective advantage of \nless-loaded inversions", trans = "log10", limits=c(0.000009, max(ResultHeatNonNeutral$CumulSelectAdv)), breaks=c(0.01,0.001, 0.0001, 0.00001))+
  ThemeSobr+
  theme(panel.border = element_blank(),
        # legend.position = c(0.5,1.1),
        # legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        strip.text = element_text(face="bold"),
        legend.background = element_blank(),
        plot.title=element_text(size=11, face="bold",hjust=0.5),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        plot.margin = margin(29, 3, 3, 5, "pt"))#+
SelectAdv

base=ggplot(ResultHeatNonNeutral[ResultHeatNonNeutral$U %in% c(0.05, 0.01, 0.001),])
# Panel A
LessLoadProb=base+geom_tile(aes(x=h, y=s, fill=TotProb))+
  ylab(expression(paste("Selection coefficient (",italic("s"), ")")))+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))+
  # coord_cartesian(ylim=c(0,0.1), expand = F)+
  facet_grid(. ~ paste0("U=", U))+
  scale_fill_viridis(option = "B", "Probability of being less loaded\nfor an inversion")+
  ThemeSobr+
  theme(panel.border = element_blank(),

        legend.text = element_text(face="bold"),
        strip.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        legend.background = element_blank(),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        plot.title=element_text(size=11, face="bold",hjust=0.5),
        )

LessLoadProb

PlotDeter=plot_grid(LessLoadProb, SelectAdv, nrow=2, labels=c('a','b'), label_size = 25)
PlotDeter


## Figure 2, panel C ##
#### Reanalyses of the dataset from Jay et al. 2022 ### 
### The data can be downloaded here : https://doi.org/10.6084/m9.figshare.19961033
Simul=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/InversionTrajectories_N=1000_Fig3.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000

#Set the column names
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "StartInv", "EndInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")

# Inversions between position 1-10Mb are on sex chromosome, Inversions between position 10-20Mb are on Autosome. Modify the "Chromosome" column to indicate this
Simul[Simul$StartInv>10000000,]$Chromosome="Autosome"

# Create a new column indicating the colums size
Simul$InvSize=Simul$EndInv - Simul$StartInv 


Simul$Code=paste(Simul$N,Simul$u,Simul$r,Simul$h,Simul$s,Simul$InvSize,Simul$Chromosome, Simul$Rep, sep="_") # Define a code identifying eavh simulation
summarySub=Simul %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)
summarySub2=Simul %>% group_by(N, u, r, h, s, InvSize, Chromosome, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)

# Some stast cited in the paper
SubsetForPaperStat=summarySub2[(summarySub2$u==1e-9 & summarySub2$h==0.3 & summarySub2$s==-0.001 & summarySub2$N==1000 & summarySub2$InvSize==2000000),]
SubsetForPaperStatFixed=SubsetForPaperStat[(SubsetForPaperStat$Chromosome=="Autosome" & SubsetForPaperStat$maxGen==24991),]
SubsetForPaperStatFixed=SubsetForPaperStat[(SubsetForPaperStat$Chromosome=="Y" & SubsetForPaperStat$maxGen==24991),]
sum(SubsetForPaperStat[SubsetForPaperStat$Chromosome=="Autosome" ,]$maxGen>15500)
sum(SubsetForPaperStat[SubsetForPaperStat$Chromosome=="Autosome" ,]$maxFreq>0.95)
summarySub2Auto=summarySub2[summarySub2$Chromosome=="Autosome",]
summarySub2AutoFixed=summarySub2Auto[summarySub2Auto$maxFreq>0.95,]
sum(summarySub2AutoFixed$InitMutNumb==0)
summarySub2Y=summarySub2[summarySub2$Chromosome=="Y",]
summarySub2YFixed=summarySub2Y[summarySub2Y$maxFreq==0.25,]
sum(summarySub2YFixed$InitMutNumb==0)

FixedSimul=summarySub[summarySub$maxFreq>0.95,]$Code #Grep the code of the inversion that have fixed
LineToAdd2=Simul[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate ending states of fixed inversion.

Counter=0
for (i in FixedSimul){ 
  Counter=Counter+1
  print(Counter)
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimul=Simul[(Simul$Code==i & Simul$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimul$Gen=24991  #Change it to 24991 (last recorded generation)
  FalseEndGoodSimul$Freq=1 #Set its frequency to 1
  FalseEndGoodSimul$MeanMutInv=Simul[(Simul$Code==i & Simul$Gen==LastGen),]$MeanMutInv #Consider that its mutation number hadn't changed
  FalseEndGoodSimul$MinMutInv=Simul[(Simul$Code==i & Simul$Gen==LastGen),]$MinMutInv
  FalseMiddleGoodSimul=FalseEndGoodSimul  #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimul$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimul, FalseMiddleGoodSimul) 
}

GoodSimulComplete=rbind(Simul,LineToAdd2) #Add these false simulation end to the simulation data.frame
GoodSimulComplete[GoodSimulComplete$Chromosome=="Y",]$Freq=4*GoodSimulComplete[GoodSimulComplete$Chromosome=="Y",]$Freq # Multiply the frequency of Y inversion by 4 to get the frequency of inversion on the Y chromosome
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color Palette
GoodSimulComplete[GoodSimulComplete$Chromosome=="Autosome",]$Rep=100000+GoodSimulComplete[GoodSimulComplete$Chromosome=="Autosome",]$Rep
GoodSimulComplete$Parameter=paste0("s=", GoodSimulComplete$s, ", h=", GoodSimulComplete$h, ", u=", GoodSimulComplete$u)
options(scipen=0)
GoodSimulComplete$Gen=GoodSimulComplete$Gen-15000 #Useless here, but start at generation 0, not 15000


DataAllEnd=GoodSimulComplete %>% group_by(u, h, s, Rep, Chromosome, InvSize) %>% summarise(MaxGen=max(Gen), MaxFreq=max(Freq), InitMutInv=MinMutInv[Gen==1],  InitMutNoInv=MeanMutNoInv[Gen==1])#Grep the end generation of each simulation
DataAllEnd$state="Lost"
DataAllEnd[(DataAllEnd$MaxGen==9991 & DataAllEnd$MaxFreq==1),]$state="Fixed"
DataAllEnd[(DataAllEnd$MaxGen==9991 & DataAllEnd$MaxFreq<1),]$state="Segregating"
DataAllEnd$Parameter=paste0("s=", DataAllEnd$s, ", h=", DataAllEnd$h, ", u=", DataAllEnd$u)
write.table(DataAllEnd, "~/Project/MutationShelteringV2/Output/DataAllEnd.tsv",sep="\t", quote=F, row.names=F ) # Save dataset, to save time, just in case

DataAllEnd=read.table("~/Project/MutationShelteringV2/Output/DataAllEnd.tsv",sep="\t",header=T )
Col=scales::viridis_pal(begin=0.3, end=0.9, option="C")(5)
options(scipen = 999)
DataAllEndYFixed=DataAllEnd[DataAllEnd$Chromosome=="Y" & DataAllEnd$state=="Fixed",] # Only consider inversion fixed on the Y chromosome
DataAllEndYFixed$DiffMutNumb=DataAllEndYFixed$InitMutInv - DataAllEndYFixed$InitMutNoInv #Absolute difference in mutation load between the inversion and other haplotype, at generation 0
DataAllEndYFixed$RelatDiffMutNumb=DataAllEndYFixed$InitMutInv/DataAllEndYFixed$InitMutNoInv #Relative difference in mutation load between the inversion and other haplotype, at generation 0


DataAllEndYFixed$U=DataAllEndYFixed$InvSize * DataAllEndYFixed$u # To express result in function of U

base=ggplot(DataAllEndYFixed[DataAllEndYFixed$s<0 & DataAllEndYFixed$U %in% c(0.01, 0.001, 0.05) & DataAllEndYFixed$h>0,])
PlotSimul=base+
  geom_hline(aes(yintercept=1), linetype="dashed", color="grey")+
  geom_boxplot(aes(fill=as.factor(s), y=RelatDiffMutNumb, x=as.factor(h) ), outliers = F)+
  facet_grid(.~ paste0("U=", U))+
  scale_fill_manual(expression(paste("Selection coefficient (",italic("s"), ")")), values=Col)+
  ThemeSobr+
   theme(panel.border = element_blank(),
        #legend.position = c(0.5,1.08),
        #legend.direction = "horizontal",
        panel.grid.major = element_line(colour = "grey95"),
        axis.title.y = element_text(face="bold", size=11),
        strip.text = element_text(face="bold"),
        legend.background = element_blank(),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        plot.title=element_text(size=11, face="bold",hjust=0.5),)+
  ylab("Relative number of mutations in\nfixed inversions compared to average inversions")+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))
  #labs(title="Relative number of mutations in fixed inversions")
PlotSimul

PlotAll=plot_grid(LessLoadProb, SelectAdv, PlotSimul, nrow=3, labels=c('a','b', 'c'), label_size = 16, label_y = 1.02, label_x=-0.005)
PlotAll

save_plot("~/Paper/ModelSexChrom_V2/Plots/Figure2_24-02-2025.png", PlotAll, ncol=2, nrow=3)
save_plot("~/Paper/ModelSexChrom_V2/Plots/Figure2_24-02-2025.pdf", PlotAll, ncol=2, nrow=3)


### Figure 3 ###
### Panels A & D ### Deterministic simulations
library(grid)
library(gridExtra)
library(ggpubr)
### Deterministic simulation ###
library(data.table)
Line=paste("time","h","s","u","n","r","Recomb","P",
           "FXN","FXI","FXNm","FXIm","FXNf","FXIf","FXm","FXf",
           "FYN","FYI","FY","Wm","Wf","D","q", sep=" ")
u=1e-08
n=2000000
File=paste0("~/Paper/ModelSexChrom/V3/CleanDataset/TimeSimul_h_s_u=",u, "n=", n, "_XYsyst_NoMutAccumul_Example.txt")
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

### Panel A ###
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color palette
base=ggplot(Table)
PlotTrajDeter=base+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  geom_line(aes(x=time, y=FYI, linetype=as.factor(s), color=Linkage), size=1, alpha=1.0)+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_linetype_manual("s",values=c("solid", "dashed", "dotted"))+
  scale_color_manual("Linkage", values=Col, guide=F)+
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
  #geom_dl(aes(x=time, y=FYI, label=Linkage, color=Linkage),method=list('last.bumpup', cex =0.8, hjust = 0.00, vjust=0.1))
PlotTrajDeter

## Panel D ###
ResultHeatNonNeutral=read.table("~/Paper/ModelSexChrom_V2/Plots/DeterministicSimul.tsv", header=T, stringsAsFactors = F)
ResultHeatNonNeutralLong=ResultHeatNonNeutral %>% pivot_longer(cols=c("FracFixedY", "FracFixedAuto"), values_to = "FracFixed", names_to = "Chrom")

ResultHeatNonNeutralLong[ResultHeatNonNeutralLong$Chrom=="FracFixedY",]$Chrom="Y"
ResultHeatNonNeutralLong[ResultHeatNonNeutralLong$Chrom=="FracFixedAuto",]$Chrom="Autosome"

ResultHeatNonNeutralLong$s=-ResultHeatNonNeutralLong$s
base=ggplot(ResultHeatNonNeutralLong[ResultHeatNonNeutralLong$U %in% c(0.001, 0.01, 0.05),])
# Panel A
options(scipen = 999)
FracFixed=base+geom_tile(aes(x=h, y=s, fill=FracFixed))+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  # coord_cartesian(ylim=c(0,0.1), expand = F)+
  facet_grid(paste0("U=",U) ~ Chrom)+
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
  )+ggtitle("Fraction of inversions expected to fix without drift")
FracFixed


### Panel B,C, E ###
#### Reanalyses of the dataset from Jay et al. 2022 ### 
### The data can be downloaded here : https://doi.org/10.6084/m9.figshare.19961033 # Not that the data for N=10000 and U=0.001 is absent in this repository. The data will be uploaded upon publication

#Panel E
Simul10k=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/InversionTrajectories_N=10000_Fig3.txt",sep=""), stringsAsFactors = F) # Data for N=10000
colnames(Simul10k)=c("N", "u", "r", "h", "s", "Gen", "StartInv", "EndInv", "Rep", 
                     "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                     "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                     "InvFit", "NoInvFit","Freq","Chromosome")

Simul1k=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/InversionTrajectories_N=1000_Fig3.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
#Set the column names
colnames(Simul1k)=c("N", "u", "r", "h", "s", "Gen", "StartInv", "EndInv", "Rep", 
                    "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                    "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                    "InvFit", "NoInvFit","Freq","Chromosome")

Simul=full_join(Simul10k, Simul1k) # Merge the datasets 
Simul$Position="Y" #Define position 
Simul[Simul$StartInv>10000000,]$Position="Autosome" # Inversion at position > 10000000 are on the autosome
Simul[Simul$Position=="Y",]$Freq=Simul[Simul$Position=="Y",]$Freq * 4 # Define Y inversion frequency
Simul$InvSize=Simul$EndInv - Simul$StartInv #Inversion size


Summary=Simul %>% group_by(N,u,r,h,s,InvSize,Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen))
Summary$State="LostEarly"
Summary[Summary$maxGen==24991,]$State="Segregating"
Summary[Summary$maxFreq>0.95,]$State="Fixed"
Summary$StateCode=1
Summary[Summary$State=="LostEarly",]$StateCode=0
Summary[Summary$State=="Segregating",]$StateCode=0
DataSummary= Summary %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), n=n())

options(scipen=0)
DataSummary$U=DataSummary$u * DataSummary$InvSize
Col3=scales::viridis_pal(begin=0.1, end=0.8, option="A")(3)
Base1k=ggplot(DataSummary[(DataSummary$s %in% c(-0.001, -0.01, -0.1) & DataSummary$U %in% c(0.001, 0.01, 0.05) & DataSummary$N==1000),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv1k=Base1k+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  geom_point(aes(x=h,color=as.factor(s)), size=5)+
  scale_color_manual("s", values=Col3)+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))+
  ThemeSobr+
  theme(text = element_text(size=18),
        panel.grid.major = element_line(color="grey95"),
        panel.background = element_blank(),
        legend.title = element_text(size = 20),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        legend.text = element_text(size = 18),
        legend.position = "none",
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face="bold"),
        plot.title= element_text(hjust = 0.5, face="bold", size=16),
        plot.margin = unit(c(0, -3, 0, 0), "cm"))+
  facet_grid(paste0("U=", U) ~ Position, scale="free")+
  ggtitle("N=1000")
PlotPropInvSpread_AllInv1k

Base10k=ggplot(DataSummary[(DataSummary$s %in% c(-0.001, -0.01, -0.1) & DataSummary$U %in% c(0.001, 0.01, 0.05) & DataSummary$N==10000),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv10k=Base10k+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  geom_point(aes(x=h,color=as.factor(s)), size=5)+
  scale_color_manual("s", values=Col3)+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))+
  ThemeSobr+
  theme(text = element_text(size=18),
        panel.grid.major = element_line(color="grey95"),
        panel.background = element_blank(),
        legend.title = element_text(size = 20),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        legend.text = element_text(size = 18),
        plot.title= element_text(hjust = 0.5, face="bold", size=16),
        axis.title.y = element_blank(),
        strip.text.x = element_text(face="bold"),
        plot.margin = unit(c(0, 0, 0, -1), "cm"))+
  facet_grid(paste0("U=", U) ~ Position, scale="free")+
  ggtitle("N=10,000")
PlotPropInvSpread_AllInv10k

legend <- get_legend(PlotPropInvSpread_AllInv10k)

### Stats indicated in the paper
DataSummary$hCat="Inferior to 0.20"
DataSummary[DataSummary$h>=0.20,]$hCat="Superior to 0.20"

SimulSumSumWide_FracFixed=DataSummary %>%
  pivot_wider(id_cols=-c( n, hCat), values_from = ProbSpread, names_from = Position) %>%
  mutate(RatioYvsAuto=Y/Autosome) #Reshape the dataframe to compare the rate of fixation on sex chromosomes and autosomes. The column "Y" contains de rate of inversion fixation on Y chromosomes, while "Autosome" contains de rate of inversion fixation on autosomes. The Column "RatioYvsAuto" contains the ratio of these value (>1 if the rate is higher on sex chromosomes, <1 otherwise)
SimulSumSumWide_FracFixed[(SimulSumSumWide_FracFixed$Y==0 & SimulSumSumWide_FracFixed$Autosome==0),]$RatioYvsAuto=1 
SimulSumSumWide_FracFixed[(SimulSumSumWide_FracFixed$Y>0 & SimulSumSumWide_FracFixed$Autosome==0),]$RatioYvsAuto=Inf


sum(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==1000,]$RatioYvsAuto>1)/nrow(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==1000,]) #Inversion are more likely to fix on Y chromosomes than on Autosomes in 94.9 % of the parameter space when N=1000.
sum(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==10000,]$RatioYvsAuto>1)/nrow(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==10000,]) #Inversion are more likely to fix on Y chromosomes than on Autosomes in 65.8 % of the parameter space when N=10000.
sum(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==1000,]$Y)/sum(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==1000,]$Autosome) #Inversion were 6.79 more likely to fix on Y chromosomes than on Autosomes in the parameter space studied.
sum(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==10000,]$Y)/sum(SimulSumSumWide_FracFixed[SimulSumSumWide_FracFixed$N==10000,]$Autosome) #Inversion were 341 more likely to fix on Y chromosomes than on Autosomes in the parameter space studied.
sum(SimulSumSumWide_FracFixed[(SimulSumSumWide_FracFixed$h<0.20 & SimulSumSumWide_FracFixed$N==1000),]$Y)/sum(SimulSumSumWide_FracFixed[(SimulSumSumWide_FracFixed$h<0.20 & SimulSumSumWide_FracFixed$N==1000),]$Autosome) #Inversion were 14.54 more likely to fix on Y chromosomes than on Autosomes in the parameter space with relatively recessive mutations 
sum(SimulSumSumWide_FracFixed[(SimulSumSumWide_FracFixed$h<0.20 & SimulSumSumWide_FracFixed$N==10000),]$Y)/sum(SimulSumSumWide_FracFixed[(SimulSumSumWide_FracFixed$h<0.20 & SimulSumSumWide_FracFixed$N==10000),]$Autosome) #Inversion were 2295 more likely to fix on Y chromosomes than on Autosomes in the parameter space with relatively recessive mutations


#Estimate fixation time #Not used in the paper
SimulYFix=Simul[(Simul$Chromosome=="Y" & Simul$Freq==0.25),]
SimulYFixFirstGenFix=SimulYFix %>% group_by(u,r,h,s,StartInv, InvSize, Chromosome, Rep) %>% summarise(FixTime=min(Gen))
SimulAutoFix=Simul[(Simul$Chromosome=="Autosome" & Simul$Freq>0.95),]
SimulAutoFixFirstGenFix=SimulAutoFix %>% group_by(u,r,h,s,StartInv, InvSize, Chromosome, Rep) %>% summarise(FixTime=min(Gen))
Merged=rbind(SimulAutoFixFirstGenFix, SimulYFixFirstGenFix)
Merged %>% group_by(Chromosome) %>% summarise(MeanFixTime=mean(FixTime)-15000) # Difference in fixation time 

## Figure 3B-C ##
library(tidyverse)
library(cowplot)
library(ggrastr)
## Load the file containing the simulation results (The same as in Jay et al. (2022), to be downloaded here: dx.doi.org/10.6084/m9.figshare.19961033)

Simul1=Simul[(Simul$u==1e-08 & Simul$h==0.2 & Simul$s==-0.001 & Simul$InvSize==5000000 & Simul$N==1000),]
Simul2=Simul[(Simul$u==1e-08 & Simul$h==0.1 & Simul$s==-0.01 & Simul$InvSize==1000000 & Simul$N==1000),]
Simul3=Simul[(Simul$u==1e-09 & Simul$h==0.01 & Simul$s==-0.1 & Simul$InvSize==1000000 & Simul$N==1000),]

SimulSub=rbind(Simul1, Simul2, Simul3) # Merge the data of the 3 focal simulations

SimulSub$Code=paste(SimulSub$N,SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$InvSize,SimulSub$Position, SimulSub$Rep, sep="_") # Define a code identifying eavh simulation
summarySub=SimulSub %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)
FixedSimul=summarySub[summarySub$maxFreq>0.95,]$Code #Grep the code of the inversion that have fixed

LineToAdd2=SimulSub[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate ending states of fixed inversion.

for (i in FixedSimul){ #Do the same thing for fixed inversions
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSub[(SimulSub$Code==i & SimulSub$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimulSub$Gen=24991  #Change it to 24991 (last recorded generation)
  FalseEndGoodSimulSub$Freq=1 #Set its frequency to 1
  FalseEndGoodSimulSub$MeanMutInv=SimulSub[(SimulSub$Code==i & SimulSub$Gen==LastGen),]$MeanMutInv #Consider that its mutation number hadn't changed
  FalseEndGoodSimulSub$MinMutInv=SimulSub[(SimulSub$Code==i & SimulSub$Gen==LastGen),]$MinMutInv
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub  #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub) 
}

GoodSimulSubComplete=rbind(SimulSub,LineToAdd2) #Add these false simulation end to the simulation data.frame
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color Palette
GoodSimulSubComplete[GoodSimulSubComplete$Position=="Autosome",]$Rep=100000+GoodSimulSubComplete[GoodSimulSubComplete$Position=="Autosome",]$Rep
GoodSimulSubComplete$U=GoodSimulSubComplete$u*GoodSimulSubComplete$InvSize
GoodSimulSubComplete$Parameter=paste0("s=", GoodSimulSubComplete$s, ", h=", GoodSimulSubComplete$h, ", U=", GoodSimulSubComplete$U)
options(scipen=0)
GoodSimulSubComplete$Gen=GoodSimulSubComplete$Gen-15000
DataAllEnd=GoodSimulSubComplete %>% group_by(U, h, s, Rep, Position) %>% summarise(MaxGen=max(Gen), MaxFreq=max(Freq))#Grep the end generation of each simulation
DataAllEnd$state="Lost"
DataAllEnd[(DataAllEnd$MaxGen>=9991 & DataAllEnd$MaxFreq==1),]$state="Fixed"
DataAllEnd[(DataAllEnd$MaxGen>=9991 & DataAllEnd$MaxFreq<1),]$state="Segregating"
SumEnd=DataAllEnd %>% group_by(U, h, s, Position) %>% count(state, sort = TRUE)  #Summarize the data for plotting summary number
SumEnd=SumEnd[SumEnd$state!="Lost",]
SumEnd$Pos=0.10
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Position=="Y"),]$Pos=1.0 #Define Position for plotting the text
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Position=="Autosome"),]$Pos=0.9
SumEnd$Parameter=paste0("s=", SumEnd$s, ", h=", SumEnd$h, ", U=", SumEnd$U)

#Plot the data 
base=ggplot(GoodSimulSubComplete) 
### Inversion frequency
### Figure 3B####
PlotTrajDrift=base+rasterize(geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Position)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=1, size=0.2)+
  geom_text(data=SumEnd, aes(x=8500, y=Pos, label=paste0("n=",n), color=as.factor(Position)), 
            vjust = -0.5, hjust = 0, size=4, show.legend = FALSE)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-linked"))+
  xlab("")+ylab("Inversion frequency")+
  ggtitle("Trajectory of inversions with drift")+
  ThemeSobr+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=18, face="bold",hjust=0.5, vjust=2),
    strip.placement = "outside",
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face="bold"),
    legend.title = element_text(size = 18, face="bold"),
    legend.text = element_text(size = 16),
    strip.text.x = element_text(size = 14, face="bold"),
    strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
    text = element_text(size=16),
    legend.position="none",
    axis.title = element_text(size = 16, face="bold"))+
  #guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid( ~ Parameter)
PlotTrajDrift

#Mutation number
### Figure 3C####
PlotMutAccul=base+rasterize(geom_line(data=GoodSimulSubComplete[GoodSimulSubComplete$MeanMutInv>0,], aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Position)), size=0.2, alpha=0.3), dpi=300)+ #Mutation number 
  #geom_hline(yintercept = 0, linetype=1, size=0.2)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-linked"))+
  xlab("Generation")+ylab("Number of mutations\nin inversions")+
  ThemeSobr+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
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
    axis.title.y = element_text(size = 18, face="bold"),
    axis.title.x = element_text(size = 18)
  )+
  guides(color = guide_legend(override.aes = list( alpha=1)))+
  facet_wrap(. ~ Parameter, ncol = 4, scales="free_y")
PlotMutAccul

legend2=get_legend(PlotMutAccul)
PlotMutAccul=PlotMutAccul+theme(legend.position="none")


# Merge plots for making Figure 3
PlotsTrajDrift <- align_plots(PlotTrajDrift,PlotMutAccul, align = 'hv', axis = 'rltp') #Align plots B and C
MergedPlotTrajDrift=plot_grid(PlotsTrajDrift[[1]],legend2,PlotsTrajDrift[[2]], ncol=2, labels = c('b', 'c'), rel_widths = c(3,0.5,3), label_size = 25)
MergedPlotTraj=plot_grid(PlotTrajDeter,MergedPlotTrajDrift, ncol=2, labels = c('a', ''), rel_widths = c(1.5,3), label_size = 25) #Align plots A, B and C


legendDeter= get_legend(FracFixed) 
## Align plot D, E1 and E2
MergeDeter_1k_10k=align_plots(FracFixed + theme(legend.position = "none", axis.text.x = element_text(size=12)) , PlotPropInvSpread_AllInv1k + theme(axis.text.x = element_text(size=12), plot.background = element_blank()), PlotPropInvSpread_AllInv10k + theme(legend.position = "none",axis.text.x = element_text(size=12), plot.background = element_blank()) , align = 'hv', axis = 'rltp')
PlotMerged=plot_grid(MergeDeter_1k_10k[[1]],legendDeter, MergeDeter_1k_10k[[2]], MergeDeter_1k_10k[[3]], legend, rel_widths = c(1,0.1,1,1,0.1), nrow=1, labels = c('d', '', 'e', '', ''))
PlotMerged

## Title for plot E
title <- ggdraw() + 
  draw_label(
    "Fraction of inversions fixed after 10,000 generations",
    fontface = 'bold',
    size=20,
    x = 0,
    y=0,
    hjust = -1,
    vjust = 1
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, -7, 0)
  )

# Merge Plot title with plots
MergedWtTitle=plot_grid(
  title, PlotMerged,
  ncol = 1,
  rel_heights = c(0.1, 1)
)
MergedWtTitle


MergedPlot=plot_grid(MergedPlotTraj,MergedWtTitle, nrow = 2, labels = c('', ''), label_size = 25)

# Save to file #
save_plot("~/Paper/ModelSexChrom_V2/Plots/Figure3_All_V2.png", MergedPlot, ncol=2, nrow=5, base_aspect_ratio = 2)
save_plot("~/Paper/ModelSexChrom_V2/Plots/Figure3_All_V2.pdf", MergedPlot, ncol=2, nrow=5, base_aspect_ratio = 2)


### Figure 4 ###
#Analyses of new dataset 
library(cowplot)
library(ggplot2)
library(tidyverse)
library(viridis)
Col2=scales::viridis_pal(begin=0.3, end=0.8, option="C")(2)
Col3=scales::viridis_pal(begin=0.2, end=0.8, option="C")(3)
Col5=scales::viridis_pal(begin=0.1, end=0.8, option="A")(5)
themeSimple=    theme(
  panel.grid.major.y = element_line(colour = "grey", size=0.1),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text=element_text(size=16),
  axis.line = element_line(colour = "grey"),
  axis.title = element_text(face="bold"))
options(scipen=999)


df=read.table("~/Project/MutationShelteringV2/Output/IntroInv_All_Summarise.txt", stringsAsFactors = F, header = T, sep=",")
df$PropFix=df$Fixed/(df$NotFixed+df$Fixed) # Fixation rate
df$Status=paste0(df$Chrom, "\n(N=",df$N,")") 

## Creating new column for plotting labels"
df$sInvClean="Intrinsically neutral inversions"
df[df$sInv==0.01,]$sInvClean="Intrinsically beneficial inversions\n(sInv=0.01)"
df[df$sInv==0.05,]$sInvClean="Intrinsically beneficial inversions\n(sInv=0.05)"
df$U=df$mu*df$Size
df$SizeLegend="2 Mb"
df[df$Size==5000000,]$SizeLegend="5 Mb"
df$hLegend="With fully recessive mutations"
df[df$h==0.3,]$hLegend="Without fully recessive mutations"

# Changing factor order
df$sInvClean=factor(df$sInvClean, levels =  c("Intrinsically neutral inversions","Intrinsically beneficial inversions\n(sInv=0.01)","Intrinsically beneficial inversions\n(sInv=0.05)"))

base=ggplot(df[(df$sInv %in% c(0.0,0.01)),],aes(x=as.factor(s), group=as.factor(s), fill=as.factor(hLegend), y=PropFix, shape=as.factor(U), color=as.factor(hLegend)))
PlotProbFix2=base+geom_point(size=8, position=position_dodge(width = 0.50)) + 
  scale_fill_manual("Distribution of dominance coefficients", values=Col5[4:5])+
  scale_color_manual("Distribution of dominance coefficients", values=c("grey10", "grey10"))+
  scale_shape_manual("U", values = c(21,24))+
  # scale_color_manual("Selection coefficient\nof mutations (s)", values=c("grey10", "grey10"), guide=F)+
  expand_limits(y = 0)+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  themeSimple+theme(#legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    legend.title = element_text(face="bold"),
    strip.text = element_text(face="bold"),
    strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  xlab("Selection coefficient of mutations (s)")+
  facet_grid(sInvClean ~ Status, scale="free")
PlotProbFix2

save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots/Figure4.svg"), PlotProbFix2, nrow=2, ncol=2)    

dfSub=df[df$sInv<0.05,]
dfSum=dfSub %>% group_by(Status, sInv) %>% summarise( n=sum(n), Fixed=sum(Fixed), NotFixed=sum(NotFixed)) # To compare fixation probility


################# Figure 6 ################
# To get the fraction of the sex chromosome not recombining,
# in the directory containing the result of all simulations, do (in a bash terminal):
# for i in *NRecomb_IndivSimulation_OnlyXY_Optim.txt ; do base=${i%%.txt} ; ../../Util/ParseRecombinationOutput_SepSex.pl -i $i -o $base.parsed.txt ; done
# for i in *_NRecomb_IndivSimulation_OnlyXY_Optim.parsed.txt ; do N_BP=`echo $i | sed 's/.*-BP=//' | sed 's/_.*//'`; echo $N_BP ;  REP=`echo $i | sed 's/.*_Rep_//' | sed 's/_.*//'`; echo $REP; echo -ne "$N_BP\t$REP\t" >> LastGenRecomb.txt ; tail -n 1 $i >> LastGenRecomb.txt ; done
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/LastGenRecomb_Fig6.txt", stringsAsFactors = F, header=F)
Pos=seq(1:(length(colnames(Data))-3))
colnames(Data)=c("N_BP","Rep","Generation", Pos) # Set the position along the genome in Mb
Data$NumberNonRecomb=0 #Number of Mb with no recombination
for (i in 1:nrow(Data)){ #For all simulation
  for (o in 4:103){ #For all position on the sex-chromosome
    if (Data[i,o]==0.0){ #If no recombination is recorded
      Data$NumberNonRecomb[i]=Data$NumberNonRecomb[i]+1 #Increment the number of non-recombining region by 1
    }
  }
}

Data$FracNonRecom=Data$NumberNonRecomb/100 #Define here the fraction of the sex-chromosme not recombining (the sex chromosome is 100Mb long)
Col4=scales::viridis_pal(begin=0.3, end=0.8, option="B")(4) #Color Palette
# Panel A #
baseNRecomb=ggplot(Data)
PlotNRecomb=baseNRecomb + geom_boxplot(aes(x=as.factor(N_BP), y=FracNonRecom, fill=as.factor(N_BP)), outlier.shape=NA) +
  geom_jitter(aes(x=as.factor(N_BP), y=FracNonRecom, fill=as.factor(N_BP)), size=2, alpha=0.4, shape=21, width = 0.1)+
  scale_fill_manual(values = Col4, guide=F)+
  ThemeSobr+
  theme(
        plot.margin = margin(3, 3, 3, 3, "pt"))+
  ylab("Fraction of the sex chromosome \n not recombining")

#Panel B#
# To get a summary of the inversion and reversion that have appeared, 
# in the directory containing the result of all simulations, do (in a bash terminal):
# for i in *_NewInv_Optim.txt ; do base=${i%%_NewInv_Optim.txt}; NInv=` cat $i | wc -l`; RevFile=$base"_Reversion_IndivSimulation_OnlyXY_NbMut_Optim.txt"; if [ -f "$RevFile" ]; then NRev=`cat $RevFile | wc -l` ; else NRev=0; fi; BlockRevFile=$base"_BlockedReversion_IndivSimulation_OnlyXY_NbMut_Optim.txt"; if [ -f "$BlockRevFile" ]; then NBlockRev=`cat $BlockRevFile | wc -l` ; else NBlockRev=0; fi;  N_BP=`echo $i | sed 's/.*-BP=//' | sed 's/_.*//'`; echo $N_BP ;  REP=`echo $i | sed 's/.*_Rep_//' | sed 's/_.*//'`; echo -ne "$N_BP\t$REP\t$NRev\t$NInv\t$NBlockRev\n" >> Nreversion_N=1000.txt ; done
DataRev=read.table("~/Paper/ModelSexChrom_V2/Datasets/Nreversion_N=1000_Fig6.txt", stringsAsFactors = F, header=F)
colnames(DataRev)=c("N_BP","Rep","N_Revers", "N_Inv", "NBlockRev")
baseNRev=ggplot(DataRev)
PlotNRev=baseNRev + geom_boxplot(aes(x=as.factor(N_BP), y=N_Revers, fill=as.factor(N_BP)), outlier.shape=NA) +
  geom_jitter(aes(x=as.factor(N_BP), y=N_Revers, fill=as.factor(N_BP)), size=2, alpha=0.4, shape=21, width = 0.1)+
  scale_fill_manual(values = Col4, guide=F)+
  scale_y_log10()+
  ThemeSobr+xlab("Number of potential inversion breakpoints (k)")+
  theme(plot.margin = margin(3, 30, 3, 3, "pt"))+
  ylab("Number of reversions (log scale)")
PlotNRev

plotsAB <- align_plots(PlotNRecomb, PlotNRev,align = 'hv', axis = 'rltp') #Align plots
MergedPlot=plot_grid(plotsAB[[1]], plotsAB[[2]], ncol=2, labels = c('a', 'b'))
save_plot("~/Paper/ModelSexChrom_V2/Plots/Fig6.png", MergedPlot, nrow=1, ncol=2, base_aspect_ratio = 1.1)
save_plot("~/Paper/ModelSexChrom_V2/Plots/Fig6.pdf", MergedPlot, nrow=1,ncol=2, base_aspect_ratio = 1)

