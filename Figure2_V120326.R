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
s_vals <- (seq(sqrt(0.001), sqrt(0.1), length.out = 100))^2
for (U in c(0.001, 0.01, 0.05))
{
  n=as.integer(U/u)
  for (h in seq(0.001,0.5,0.0025))
  {
    for (s in s_vals) 
    {
      #print(paste(c(U, h, s)))
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

write.table(ResultHeatNonNeutral,"~/Project/MutationShelteringV2/Revision1/Output/DeterministicSimul_sqrt.tsv", sep="\t", row.names = F, quote = F)
ResultHeatNonNeutral=read.table("~/Project/MutationShelteringV2/Revision1/Output/DeterministicSimul_sqrt.tsv", header=T, stringsAsFactors = F)
mean(ResultHeatNonNeutral[ResultHeatNonNeutral$U==0.001,]$TotProb)
Subset=ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.001 ),]
mean(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.05 ),]$CumulSelectAdv)


mean(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.001 ),]$CumulSelectAdv)
mean(ResultHeatNonNeutral[(ResultHeatNonNeutral$s==0.001 ),]$CumulSelectAdv)

mean(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.05 & ResultHeatNonNeutral$s == 0.001 & (ResultHeatNonNeutral$h > 0.3 & ResultHeatNonNeutral$h < 0.31)),]$CumulSelectAdv)
quantile(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.05 & ResultHeatNonNeutral$s == 0.010123967 & (ResultHeatNonNeutral$h > 0.3 & ResultHeatNonNeutral$h < 0.32)),]$CumulSelectAdv)

sum(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.05),]$CumulSelectAdv > 0.001)
sum(ResultHeatNonNeutral[(ResultHeatNonNeutral$U==0.05),]$CumulSelectAdv > 0.00000001)

ResultHeatNonNeutral$s=-ResultHeatNonNeutral$s
options(scipen = 9)


signed_sqrt_trans <- function() {
  trans_new(
    "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse   = function(x) sign(x) * (x^2),
    domain = c(-Inf, Inf)
  )
}

base=ggplot(ResultHeatNonNeutral[ResultHeatNonNeutral$U %in% c(0.05, 0.01, 0.001),])
# Panel B
SelectAdv=base+geom_tile(aes(x=h, y=s, fill=CumulSelectAdv))+
  ylab(expression(paste("Selection coefficient (",italic("s"), ")")))+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))+
  # coord_cartesian(ylim=c(0,0.1), expand = F)+
  facet_grid(
    . ~ paste0("U[L]==", U),
    labeller = label_parsed
  )+
  scale_fill_viridis(option = "B", "Average selective advantage of \nless-loaded inversions", trans = "log10", limits=c(0.000009, max(ResultHeatNonNeutral$CumulSelectAdv)), breaks=c(0.01,0.001, 0.0001, 0.00001))+
  scale_y_continuous(trans = signed_sqrt_trans(), breaks=c(-0.001, -0.004, -0.01, -0.02, -0.05, -0.1))+
  ThemeSobr+
  theme(panel.border = element_blank(),
        # legend.position = c(0.5,1.1),
        # legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        #strip.text = element_text(face="bold"),
        strip.text.x = element_blank(),
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
  facet_grid(
    . ~ paste0("U[L]==", U),
    labeller = label_parsed
  )+
  scale_fill_viridis(option = "B", "Probability of being less loaded\nfor an inversion")+
  scale_y_continuous(trans = signed_sqrt_trans(), breaks=c(-0.001, -0.004, -0.01, -0.02, -0.05, -0.1))+
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

##
Simul=read.table(paste("~/Project/MutationShelteringV2/Revision1/Output/IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat_Full.24Col.txt",sep=""), stringsAsFactors = F) # Data for N=10000
colnames(Simul)=c("N", "u", "r", "h", "s", "sInv", "Gen", "StartInv", "EndInv", "Rep" ,"Chromosome",
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")
SimulCp=Simul
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
Col3=scales::viridis_pal(begin=0.20, end=0.8, option="A")(3)
options(scipen = 999)
DataAllEnd$RelatDiff=DataAllEnd$InitMutNumb/DataAllEnd$InitMutNumbNoInv

DataAllEndYFixed=DataAllEnd[(DataAllEnd$Chromosome=="Ychromosome" & DataAllEnd$State=="Fixed"),] # Only consider inversion fixed on the Y chromosome
DataAllEndYFixed$U=DataAllEndYFixed$InvSize * DataAllEndYFixed$u # To express result in function of U
DataAllEndYFixed$log2FC=log2(DataAllEndYFixed$InitMutNumb/DataAllEndYFixed$InitMutNumbNoInv)
DataAllEndYFixed$RelatDiff=DataAllEndYFixed$InitMutNumb/DataAllEndYFixed$InitMutNumbNoInv

DataAllEndYFixed$sNorm=factor(DataAllEndYFixed$sNorm, levels = c(-1,-10,-100))
base=ggplot(DataAllEndYFixed[(DataAllEndYFixed$s<0 & DataAllEndYFixed$U %in% c(0.01, 0.001, 0.05) & DataAllEndYFixed$h>0 & DataAllEndYFixed$sNorm %in% c(-1, -10, -100)),])
PlotSimul=base+
  geom_hline(aes(yintercept=1), linetype="dashed", color="grey")+
  geom_boxplot(aes(fill=sNorm, y=RelatDiff, x=as.factor(h) ), outliers = F)+
  facet_grid(
    . ~ paste0("U[L]==", U),
    labeller = label_parsed
  )+
  scale_fill_manual(expression(paste("Selection coefficient (",italic("Ns"), ")")), values=Col3)+
  ThemeSobr+
  theme(panel.border = element_blank(),
        # legend.position = c(0.5,1.1),
        # legend.direction = "horizontal",
        panel.grid.major = element_line(colour="grey97"),
        legend.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        #strip.text = element_text(face="bold", size=14),
        strip.text.x = element_blank(),
        legend.background = element_blank(),
        plot.title=element_text(size=11, face="bold",hjust=0.5),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
        axis.title.y = element_text(size=14),)+
        #plot.margin = margin(29, 3, 3, 10, "pt"))+
  ylab("Relative number of mutations in\nfixed inversions compared to average")+
  xlab(expression(paste("Dominance coefficient (",italic("h"), ")")))
#labs(title="Relative number of mutations in fixed inversions")
PlotSimul

# PlotAll=plot_grid(LessLoadProb, LessLoadProb, PlotSimul, nrow=3, labels=c('a','b', 'c'), label_size = 16, label_y = 1.02, label_x=-0.005)
PlotAll=LessLoadProb / SelectAdv / PlotSimul + plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(size = 16, face="bold"))
PlotAll
save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/Figure2_V090426.pdf"), PlotAll, nrow=3, ncol=2)    
