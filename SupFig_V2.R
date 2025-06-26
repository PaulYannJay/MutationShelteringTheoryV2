#### Papier Sheltering, V3 (2025) ## SUP FIGURE ##

library(cowplot)
library(ggplot2)
library(tidyverse)
library(viridis)
library(directlabels)
library(ggrastr)
library(ggnewscale)
library(gghalves)
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


#####################################################
################@ Figure S1 #########################
#####################################################

u=1e-08
h=0.1
s=0.01
q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
n=2000000
BinomProba=data.frame(m=double(),u=double(), Proba=double()) #Dataset 
for (m in seq(0,n*q*8)) #Only consider m varying from 0 to 8nq to save time
{
  Prob=dbinom(m,n,q) #Calculate binomial probabilities
  BinomProba[nrow(BinomProba)+1,]=c(m,u,Prob) #Fill the table (ugly way of doing that...)
}

dfProb_areaStep <- bind_rows(old = BinomProba, 
                             new = BinomProba %>% mutate(Proba = lag(Proba)),
                             .id = "source") %>% arrange(m, source)

dfProb_areaStep$Proba[1]=0
B3=n*q
B1=(q*h*n)/(1-h)
dfProb_areaStep$m=dfProb_areaStep$m-0.5 # to center values
Col=scales::viridis_pal(begin=0.2, end=0.8, option="A")(3)
base=ggplot(BinomProba, aes(x=m, y=Proba))
PlotProbaStat=base+geom_ribbon(data = dfProb_areaStep[dfProb_areaStep$m<B1+1,],aes(xmin = -1, xmax=B1, ymin = 0, ymax = Proba, fill="WII>WNI (inversion fixation)"))+
  geom_ribbon(data = dfProb_areaStep[(dfProb_areaStep$m>B1 & dfProb_areaStep$m<B3+1),],aes(xmin = -1, xmax=B3+2, ymin = 0, ymax = Proba, fill="WNI>WII (inversion maintained at intermediate frequency)"))+
  geom_ribbon(data = dfProb_areaStep[dfProb_areaStep$m>=B3,],aes(xmin = -1, xmax=max(dfProb_areaStep$m), ymin = 0, ymax = Proba, fill="WNI<WNN (inversion lost)"))+
  xlim(-0.5,n*q*3)+
  scale_fill_manual("Inversion state", values=Col)+
  annotate("segment", x=B1, xend=B1,y=0, yend=max(BinomProba$Proba)+0.05*max(BinomProba$Proba), color="black")+
  annotate("segment", x=B3, xend=B3,y=0, yend=max(BinomProba$Proba)+0.05*max(BinomProba$Proba), color="black")+
  xlab("m (number of mutations in the inversion)")+
  ylab("Probability of occurence")+
  annotate("text", x=B1, y=max(BinomProba$Proba)+0.1*max(BinomProba$Proba), label="m=qhn/(1-h)")+
  annotate("text", x=B3, y=max(BinomProba$Proba)+0.1*max(BinomProba$Proba), label="m=nq")+
  ThemeSobr+
  theme(legend.text = element_text(size=12, face="bold"),
        legend.position = c(0.75, 0.5),
        legend.title = element_text(size=12, face="bold"))+
  ggtitle(paste0("n=",n,", u=", u, ", s=", s, ", h=", h))

save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS1.png"),PlotProbaStat, nrow=1, ncol=1, base_aspect_ratio = 3) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS1.pdf"),PlotProbaStat, nrow=1, ncol=1, base_aspect_ratio = 3) #Fig. S?

#####################################################
################@ Figure S2-3 #######################
##################################################### 
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/Linked_2MbInv_Gen1000_FigS3.txt", stringsAsFactors = F) #FigS3
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/Unlinked_2MbInv_Gen1000_FigS2.txt", stringsAsFactors = F) #FigS2
colnames(Data)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Data=Data[Data$MutInv != -1, ]
Data$Code=paste(Data$u,Data$r,Data$h,Data$s,Data$DebInv,Data$FinInv, Data$Rep, sep="_")
Data$MutInv=Data$MutInv ### To remove from count the inversion mutation
DataBegin=Data[Data$Gen==15002,]
DataEnd=Data[Data$Gen==16002,]
DataBeginLost=DataBegin[! DataBegin$Code %in% DataEnd$Code, ]
FalseDataEnd=DataBeginLost
FalseDataEnd$Gen=16002
FalseDataEnd$Freq=0.0
DataAllEnd=rbind(DataEnd, FalseDataEnd)
DataAllEnd$State="Lost"
DataAllEnd[DataAllEnd$Freq>0,]$State="Segregate"
DataAllEnd$RelativeMutNumber=DataAllEnd$MutInv/DataAllEnd$MutNoInv
DataAllEnd$RelativeMutFreq=DataAllEnd$FreqMutInv/DataAllEnd$FreqMutNoInv


FocS=c(-0.01,-0.05,-0.1)
FocH=c(0.001,0.01,0.1)
DataAllEnd=DataAllEnd[(DataAllEnd$s %in% FocS & DataAllEnd$h %in% FocH),]
Col=scales::viridis_pal(begin=0.4, end=0.6, option="A")(2)
Plot=DataAllEnd %>% 
  group_split( s) %>% 
  map( ~ggplot(., aes(x=State, y=RelativeMutNumber)) +
         geom_hline(yintercept = 1.0,linetype = 2)+
         geom_half_point_panel(aes(fill=RelativeMutFreq), shape=21,size=2, range_scale = 1.0)+
         scale_fill_gradient2("Relative frequency of\nmutations in inversions",
                              low = "navyblue", 
                              mid = "white", 
                              high = "firebrick", 
                              midpoint = 1,
         ) +
         new_scale("fill") +
         geom_half_boxplot(aes(fill=State),side = "l",  errorbar.draw = FALSE, outlier.shape = NA, alpha=1)+
         scale_fill_manual("", values=Col,guide=FALSE)+
         xlab("Inversion state after 1000 generations")+
         ylab("Relative number of mutations")+
         facet_grid(paste0("s=",s)~paste0("h=", h), scales = "free_y", labeller = function(x) label_value(x, multi_line = FALSE))+
         ggtitle("Dominance")+
         theme(panel.border = element_blank(),  
               legend.background = element_blank(),
               panel.background = element_blank(),
               text = element_text(size=14),
               axis.line = element_line(colour = "black"),
               axis.title.y = element_text(face="bold"),
               axis.text.x=element_text(face="bold"),
               legend.text = element_text(size = 10),
               legend.title = element_text(size = 12),
               legend.key.size = unit(15, 'pt'),
               plot.title = element_text(hjust = 0.5, size=12, vjust=-1, face="bold")
         )
  ) %>% plot_grid(plotlist = ., align = 'hv', ncol = 1)

#save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS3.png", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)
#save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS3.pdf", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS2.png", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS2.pdf", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)

#####################################################
################@ Figure S4 #########################
#####################################################
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                      FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double(), Chrom=character())
write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_XY_Y.txt", append = F, quote=F, row.names = F)
P=0.80
Chrom="Y-linked"
for (R in c(0.0, 0.001, 0.01, 0.05, 0.5))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.1))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                            FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                            FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double())
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXNm=1.00
      FXIm=0.00
      FXNf=1.00
      FXIf=0.00
      FYN=0.99
      FYI=0.01
      FY=FYN+FYI
      FXm=FXNm+FXIm
      FXf=FXNf+FXIf
      FXI=(2/3)*FXIf + (1/3)*FXIm
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
      D=FXI*FYN - FXN*FYI
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        FYI=(TableLinkR$FYI[time-1]*TableLinkR$FXIf[time-1]*WII +
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*(1-R) + 
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FYN=(TableLinkR$FYN[time-1]*TableLinkR$FXNf[time-1]*WNN +
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*(1-R) + 
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        
        FXIm=(TableLinkR$FXIf[time-1]*TableLinkR$FYI[time-1]*WII +
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*(1-R) + 
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXNm=(TableLinkR$FXNf[time-1]*TableLinkR$FYN[time-1]*WNN +
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*(1-R) + 
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXIf=(TableLinkR$FXIf[time-1]*TableLinkR$FXIm[time-1]*WII +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXNf=(TableLinkR$FXNf[time-1]*TableLinkR$FXNm[time-1]*WNN +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXI=(2/3)*FXIf + (1/3)*FXIm
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      TableLinkR$Chrom=Chrom
      write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_XY_Y.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                      FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double(), Chrom=character())
write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_XY_X.txt", append = F, quote=F, row.names = F)
P=0.80
Chrom="X-linked"
for (R in c(0.0, 0.001, 0.01, 0.05, 0.5))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.1))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                            FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                            FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double())
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXNm=0.99
      FXIm=0.01
      FXNf=0.99
      FXIf=0.01
      FYN=1.00
      FYI=0.00
      FY=FYN+FYI
      FXm=FXNm+FXIm
      FXf=FXNf+FXIf
      FXI=(2/3)*FXIf + (1/3)*FXIm
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
      D=FXI*FYN - FXN*FYI
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        FYI=(TableLinkR$FYI[time-1]*TableLinkR$FXIf[time-1]*WII +
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*(1-R) + 
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FYN=(TableLinkR$FYN[time-1]*TableLinkR$FXNf[time-1]*WNN +
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*(1-R) + 
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        
        FXIm=(TableLinkR$FXIf[time-1]*TableLinkR$FYI[time-1]*WII +
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*(1-R) + 
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXNm=(TableLinkR$FXNf[time-1]*TableLinkR$FYN[time-1]*WNN +
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*(1-R) + 
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXIf=(TableLinkR$FXIf[time-1]*TableLinkR$FXIm[time-1]*WII +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXNf=(TableLinkR$FXNf[time-1]*TableLinkR$FXNm[time-1]*WNN +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXI=(2/3)*FXIf + (1/3)*FXIm
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      TableLinkR$Chrom=Chrom
      write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_XY_X.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}


Data_Y=read.table("~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_XY_Y.txt", stringsAsFactors = F, header = T)
Data_X=read.table("~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_XY_X.txt", stringsAsFactors = F, header = T)
FocR=c(0,0.1,1,5,50) # Value of r to focus on
FocS=c(0.01,0.1) # Value of s to focus on
Data_Y=Data_Y[Data_Y$s %in% FocS,]
Data_Y$r=Data_Y$r*100
Data_Y=Data_Y[Data_Y$r %in% FocR,]
Data_Y$Chrom="Y-linked"
Data_X=Data_X[Data_X$s %in% FocS,]
Data_X$r=Data_X$r*100
Data_X=Data_X[Data_X$r %in% FocR,]
Data_X$Chrom="X-linked"
Data=rbind(Data_X,Data_Y)
Data$Freq=Data$FYI+Data$FXI

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(5)
base=ggplot(Data[Data$Chrom=="Y-linked",])
Plot0.8_Y=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=Freq, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(),  
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.08,1.1))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))

base=ggplot(Data[Data$Chrom=="X-linked",])
Plot0.8_X=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=Freq, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(),  
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.08,1.1))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))

Plot=plot_grid(Plot0.8_Y, Plot0.8_X, nrow=2, labels=c('a','b'), label_size = 25)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS4.png", Plot, ncol = 3, nrow=4)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS4.pdf", Plot, ncol = 3, nrow=4)

#####################################################
################@ Figure S5 #########################
#####################################################

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(),
                      FYN=double(),FYI=double(), FY=double(), 
                      W=double(),D=double(), q=double())
write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_2MT.txt", append = F, quote=F, row.names = F)

for (P in c(0.8, 0.9, 0.99))
{
  for (R in c(0.0, 0.001, 0.01, 0.05, 0.5))
  {
    for (h in c(0.001, 0.01, 0.1))
    {
      for (s in c(0.01, 0.1))
      {
        TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                              FXN=double(),FXI=double(),FX=double(),
                              FYN=double(),FYI=double(), FY=double(), 
                              W=double(),D=double(),q=double())
        u=1e-08
        n=2000000
        q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
        m=floor(P*n*q)
        WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
        WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
        WII=(1-s)^m #Fitness of individual homozygous for the inversion.
        FXN=1.0
        FXI=0.00
        FYN=0.99
        FYI=0.01
        FY=FYN+FYI
        FX=FXN+FXI
        W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
        D=FXI*FYN - FXN-FYI
        FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=FirstGen
        for (time in seq(2,10000,1))
        {
          FXI2=(FXI*FYI*WII + FXI*FYN*WNI*(1-R) + FXN*FYI*WNI*R)/W
          FXN2=(FXN*FYN*WNN + FXN*FYI*WNI*(1-R) + FXI*FYN*WNI*R)/W
          FYI2=(FYI*FXI*WII + FYI*FXN*WNI*(1-R) + FYN*FXI*WNI*R)/W
          FYN2=(FYN*FXN*WNN + FYN*FXI*WNI*(1-R) + FYI*FXN*WNI*R)/W
          FXI=FXI2
          FXN=FXN2
          FYI=FYI2
          FYN=FYN2
          FY=FYN+FYI
          FX=FXN+FXI
          D=FXI*FYN - FXN-FYI
          W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
          Gen=c(time,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
          TableLinkR[nrow(TableLinkR)+1,]=Gen
        }
        write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_2MT.txt", append = T, quote=F, row.names = F, col.names = F)
      }
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(5)
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_2MT.txt", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1)
Data=Data[Data$s %in% FocS,]
Data$r=Data$r*100
FoccM=c(0, 0.1, 1, 5, 50)
Data=Data[Data$r %in% FoccM,]
base=ggplot(Data[Data$P==0.8,])
Plot0.8=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))


base=ggplot(Data[Data$P==0.9,])
Plot0.9=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

base=ggplot(Data[Data$P==0.99,])
Plot0.99=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

Plot=plot_grid(Plot0.8,Plot0.9,Plot0.99, nrow=3, labels = c('a', 'b', 'c'), label_size = 30)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS5.png",Plot, ncol=3, nrow=6)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS5.pdf",Plot, ncol=3, nrow=6)

#####################################################
################@ Figure S6 #########################
#####################################################

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(),
                      FYN=double(),FYI=double(), FY=double(), 
                      W=double(),D=double(), q=double(), Des=double(),  InitFreq=double())
write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_2MT_VarLD.txt", append = F, quote=F, row.names = F)
P=0.80
for( Des in c(1,0.5)) #When Des=1, the inversion is introduced only on Y chromosomes. When Des=0.5, the inversion is introduce in equal proportion in X and Y chromosomes
{
  for ( InitFreq in c(0.01, 0.001, 0.0001)) #Variation in the initial frequency of the inversion
  {
    for (R in c(0.0, 0.001, 0.005))
    {
      for (h in c(0.01))
      {
        for (s in c(0.1, 0.01))
        {
          TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                                FXN=double(),FXI=double(),FX=double(),
                                FYN=double(),FYI=double(), FY=double(), 
                                W=double(),D=double(),q=double(), Des=double(),  InitFreq=double())
          
          u=1e-08
          n=2000000
          q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
          m=floor(P*n*q)
          WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
          WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
          WII=(1-s)^m #Fitness of individual homozygous for the inversion.
          XRelatF=(1-Des)*InitFreq
          YRelatF=Des*InitFreq
          FXN=1.0-XRelatF
          FXI=XRelatF
          FYN=1.0-YRelatF
          FYI=YRelatF
          FY=FYN+FYI
          FX=FXN+FXI
          W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
          D=FXI*FYN - FXN*FYI
          FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q, Des, InitFreq)
          TableLinkR[nrow(TableLinkR)+1,]=FirstGen
          for (time in seq(2,10000,1))
          {
            FXI2=(FXI*FYI*WII + FXI*FYN*WNI*(1-R) + FXN*FYI*WNI*R)/W
            FXN2=(FXN*FYN*WNN + FXN*FYI*WNI*(1-R) + FXI*FYN*WNI*R)/W
            FYI2=(FYI*FXI*WII + FYI*FXN*WNI*(1-R) + FYN*FXI*WNI*R)/W
            FYN2=(FYN*FXN*WNN + FYN*FXI*WNI*(1-R) + FYI*FXN*WNI*R)/W
            FXI=FXI2
            FXN=FXN2
            FYI=FYI2
            FYN=FYN2
            FY=FYN+FYI
            FX=FXN+FXI
            D=FXI*FYN - FXN*FYI
            W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
            Gen=c(time,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q,Des,InitFreq)
            TableLinkR[nrow(TableLinkR)+1,]=Gen
          }
          write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_2MT_VarLD.txt", append = T, quote=F, row.names = F, col.names = F)
        }
      }
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(3)
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_2MT_VarLD.txt", stringsAsFactors = F, header = T)
Data$r=Data$r*100
Data=Data[Data$Des %in% c(0.5,1.0),]
Data$InitD=0
for (H in unique(Data$h))
{
  for (S in unique(Data$s))
  {
    for (I in unique(Data$InitFreq))
    {
      for (D in unique(Data$Des))
      {
        for (R in unique(Data$r))
        {
          InitDFix=Data[(Data$s==S & Data$h==H & Data$InitFreq==I & Data$Des==D & Data$time==1),]$D
          Data[(Data$s==S & Data$h==H & Data$InitFreq==I & Data$Des==D),]$InitD=InitDFix
        }
      }
    }
  }
}

base=ggplot(Data)
Plot0.8_VarFreq=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~ fct_relevel(paste0("Initial D=",-InitD), "D=0", "D=0.0001", "D=0.001", "D=0.01"))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

base=ggplot(Data)
Plot0.8_VarLD=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=-D, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~ fct_relevel(paste0("Initial D=",-InitD), "D=0", "D=0.0001", "D=0.001", "D=0.01"))+
  ylab("Linkage desequilibrium")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

Plot=plot_grid(Plot0.8_VarFreq, Plot0.8_VarLD, nrow=2, labels=c("A", "B"), label_size = 30)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS6.png", Plot, ncol=3, nrow=4)

#####################################################
################@ Figure S7 #########################
#####################################################

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(),P=double(), 
                      FXN=double(),FXI=double(),FX=double(),
                      FYN=double(),FYI=double(), FY=double(),
                      FZN=double(),FZI=double(), FZ=double(),
                      FXIYI=double(),FXNYI=double(),FXIYN=double(),FXNYN=double(),
                      FYIZI=double(),FYNZI=double(),FYIZN=double(),FYNZN=double(),
                      FXIZI=double(),FXNZI=double(),FXIZN=double(),FXNZN=double(),
                      FT=double(),W=double(),q=double())

write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_3MT", append = F, quote=F, row.names = F)
P=0.8
for (R in c(0.0, 0.001, 0.01, 0.05, 0.1, 0.5))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.1))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(),r=double(), P=double(), 
                            FXN=double(),FXI=double(),FX=double(),
                            FYN=double(),FYI=double(), FY=double(),
                            FZN=double(),FZI=double(), FZ=double(),
                            FXIYI=double(),FXNYI=double(),FXIYN=double(),FXNYN=double(),
                            FYIZI=double(),FYNZI=double(),FYIZN=double(),FYNZN=double(),
                            FXIZI=double(),FXNZI=double(),FXIZN=double(),FXNZN=double(),
                            FT=double(),W=double(),q=double())
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      F_XN=0.99
      F_XI=0.01
      F_YN=1.00
      F_YI=0.00
      F_ZN=1.00
      F_ZI=0.00
      FY=1/3
      FX=1/3
      FZ=1/3
      FXN=F_XN*FX
      FXI=F_XI*FX
      FYN=F_YN*FY
      FYI=F_YI*FY
      FZN=F_ZN*FZ
      FZI=F_ZI*FZ
      FXIYI=FXI*(FYI/(FY+FZ))+FYI*(FXI/(FX+FZ))
      FXIYN=FXI*(FYN/(FY+FZ))+FYN*(FXI/(FX+FZ))
      FXNYI=FXN*(FYI/(FY+FZ))+FYI*(FXN/(FX+FZ))
      FXNYN=FXN*(FYN/(FY+FZ))+FYN*(FXN/(FX+FZ))
      FYIZI=FZI*(FYI/(FY+FX))+FYI*(FZI/(FX+FZ))
      FYIZN=FZI*(FYN/(FY+FX))+FYN*(FZI/(FX+FZ))
      FYNZI=FZN*(FYI/(FY+FX))+FYI*(FZN/(FX+FZ))
      FYNZN=FZN*(FYN/(FY+FX))+FYN*(FZN/(FX+FZ))
      FXIZI=FZI*(FXI/(FY+FX))+FXI*(FZI/(FY+FZ))
      FXIZN=FZN*(FXI/(FY+FX))+FXI*(FZN/(FY+FZ))
      FXNZI=FZI*(FXN/(FY+FX))+FXN*(FZI/(FY+FZ))
      FXNZN=FZN*(FXN/(FY+FX))+FXN*(FZN/(FY+FZ))
      Sum=FXIYI+FXNYI+FXIYN+FXNYN+
        FYIZI+FYNZI+FYIZN+FYNZN+
        FXIZI+FXNZI+FXIZN+FXNZN
      FT=FX+FY+FZ
      W=FXIYI*WII+FXNYI*WNI+FXIYN*WNI+FXNYN*WNN+
        FYIZI*WII+FYNZI*WNI+FYIZN*WNI+FYNZN*WNN+
        FXIZI*WII+FXNZI*WNI+FXIZN*WNI+FXNZN*WNN
      
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,FZN,FZI,FZ,
                 FXIYI,FXNYI,FXIYN,FXNYN,FYIZI,FYNZI,FYIZN,FYNZN,FXIZI,FXNZI,FXIZN,FXNZN,
                 FT,W,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        
        FXI=((TableLinkR$FXIYI[time-1]+TableLinkR$FXIZI[time-1])*WII+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FXIZN[time-1])*(1-R)*WNI+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FXNZI[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FXN=((TableLinkR$FXNYN[time-1]+TableLinkR$FXNZN[time-1])*WNN+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FXNZI[time-1])*(1-R)*WNI+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FXIZN[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FYI=((TableLinkR$FXIYI[time-1]+TableLinkR$FYIZI[time-1])*WII+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FYIZN[time-1])*(1-R)*WNI+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FYNZI[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FYN=((TableLinkR$FXNYN[time-1]+TableLinkR$FYNZN[time-1])*WNN+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FYNZI[time-1])*(1-R)*WNI+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FYIZN[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FZI=((TableLinkR$FXIZI[time-1]+TableLinkR$FYIZI[time-1])*WII+
               (TableLinkR$FXNZI[time-1]+TableLinkR$FYNZI[time-1])*(1-R)*WNI+
               (TableLinkR$FXIZN[time-1]+TableLinkR$FYIZN[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FZN=((TableLinkR$FXNZN[time-1]+TableLinkR$FYNZN[time-1])*WNN+
               (TableLinkR$FXIZN[time-1]+TableLinkR$FYIZN[time-1])*(1-R)*WNI+
               (TableLinkR$FXNZI[time-1]+TableLinkR$FYNZI[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FY=FYN+FYI
        FX=FXN+FXI
        FZ=FZN+FZI
        FXIYI=FXI*(FYI/(FY+FZ))+FYI*(FXI/(FX+FZ))
        FXIYN=FXI*(FYN/(FY+FZ))+FYN*(FXI/(FX+FZ))
        FXNYI=FXN*(FYI/(FY+FZ))+FYI*(FXN/(FX+FZ))
        FXNYN=FXN*(FYN/(FY+FZ))+FYN*(FXN/(FX+FZ))
        FYIZI=FZI*(FYI/(FY+FX))+FYI*(FZI/(FX+FZ))
        FYIZN=FZI*(FYN/(FY+FX))+FYN*(FZI/(FX+FZ))
        FYNZI=FZN*(FYI/(FY+FX))+FYI*(FZN/(FX+FZ))
        FYNZN=FZN*(FYN/(FY+FX))+FYN*(FZN/(FX+FZ))
        FXIZI=FZI*(FXI/(FY+FX))+FXI*(FZI/(FY+FZ))
        FXIZN=FZN*(FXI/(FY+FX))+FXI*(FZN/(FY+FZ))
        FXNZI=FZI*(FXN/(FY+FX))+FXN*(FZI/(FY+FZ))
        FXNZN=FZN*(FXN/(FY+FX))+FXN*(FZN/(FY+FZ))
        FT=FX+FY+FZ
        W=FXIYI*WII+FXNYI*WNI+FXIYN*WNI+FXNYN*WNN+
          FYIZI*WII+FYNZI*WNI+FYIZN*WNI+FYNZN*WNN+
          FXIZI*WII+FXNZI*WNI+FXIZN*WNI+FXNZN*WNN
        
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,FZN,FZI,FZ,
              FXIYI,FXNYI,FXIYN,FXNYN,FYIZI,FYNZI,FYIZN,FYNZN,FXIZI,FXNZI,FXIZN,FXNZN,
              FT,W,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_3MT", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(6)
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_3MT", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1)
Data=Data[Data$s %in% FocS,]
FocH=c(0.001,0.01,0.1)
Data=Data[Data$h %in% FocH,]
Data$r=Data$r*100
FoccM=c(0, 0.1, 1, 5, 50)
Data=Data[Data$r %in% FoccM,]
base=ggplot(Data)
Plot0.8=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=FXI, color=paste0(r, " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(),  
        legend.position = "top",
        legend.direction = "horizontal",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=25),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.1,0.2,0.3,0.4), limits=c(0,0.5))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))#+

save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS7.png", Plot0.8, ncol = 3, nrow=2)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS7.pdf", Plot0.8, ncol = 3, nrow=2)

#####################################################
################@ Figure S8 #########################
#####################################################

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), s2=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(), FYN=double(),FYI=double(), FY=double(),
                      W=double(),D=double(),q=double())

write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_Over.txt", append = F, quote=F, row.names = F)
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
      write.table(TableLinkR, "~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_Over.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8)(6)
Data=read.table("~/Paper/ModelSexChrom_V2/Datasets/DeterministicInvEvolution_Over.txt", stringsAsFactors = F, header = T)
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

save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS8.png", Plot0.8, ncol = 3, nrow=4)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS8.pdf", Plot0.8, ncol = 3, nrow=4)


#####################################################
################@ Figure S9 #########################
#####################################################

SimulSub=read.table("~/Paper/ModelSexChrom_V2/Datasets/InversionTrajectories_N=1000_CodeB_FigureS9.txt",  header=T, stringsAsFactors = F)
SimulSubEx=subset(SimulSub, (SimulSub$h==0.1 & SimulSub$s==-0.01 & SimulSub$InvSize==2000000)) #Value to focus on
summarySub=SimulSubEx %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv))
LostSimulSub=summarySub[summarySub$maxGen<24991,]$Code # Inversion that have been lost or fixed
GoodSimulSubB=SimulSubEx #Copy, just in case
FixedSimul=summarySub[(summarySub$maxFreq>0.95 & summarySub$maxGen<24991),]$Code #Grep the code of the inversion that have fixed
LostSimul=summarySub[summarySub$maxFreq<0.95,]$Code #Mutation that have been lost
LineToAdd=SimulSubEx[0,]
for (i in LostSimul){ # As done before (Figure 3A), recreate false end state for plotting
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]
  FalseEndGoodSimulSub$Gen=24991
  FalseEndGoodSimulSub$Freq=0
  FalseEndGoodSimulSub$MeanMutInv=0
  FalseEndGoodSimulSub$MinMutInv=0
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd=rbind(LineToAdd, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}

LineToAdd2=SimulSubEx[0,]
for (i in FixedSimul){
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]
  FalseEndGoodSimulSub$Gen=24991
  FalseEndGoodSimulSub$Freq=1
  FalseEndGoodSimulSub$MeanMutInv=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]$MeanMutInv
  FalseEndGoodSimulSub$MinMutInv=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]$MinMutInv
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}


GoodSimulSubCompleteB=rbind(GoodSimulSubB,LineToAdd,LineToAdd2) #Add the false end to the dataset
GoodSimulSubCompleteB=GoodSimulSubB
GoodSimulSubCompleteB$Gen = GoodSimulSubCompleteB$Gen - 15000
options(scipen=0) #Non scientific notation
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color Palette
GoodSimulSubCompleteB$U=GoodSimulSubCompleteB$u*GoodSimulSubCompleteB$InvSize

# Panel a
base=ggplot(GoodSimulSubCompleteB[(GoodSimulSubCompleteB$Position=="Autosome"),])
FreqAuto=base+geom_line(aes(x=Gen, y=Freq, group=Rep), color=Col[1],  size=0.8, alpha=0.5)+
  facet_wrap(.~paste0("U=",U), scales = "free")+
  xlab("Generation")+
  ylab("Autosomal inversion frequency")+
  labs(title = "Inversion frequency")+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  ThemeSobr
FreqAuto

# Panel c
MinMutAuto=base+geom_line(data=GoodSimulSubCompleteB[(GoodSimulSubCompleteB$Position=="Autosome"),],
                          aes(x=Gen, y=MinMutInv, group=Rep), color=Col[1],  size=0.8, alpha=0.5)+
  facet_wrap(.~paste0("U=",U), scales = "free")+
  xlab("Generation")+
  ylab("Number of mutations")+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  labs(title = "Minimum number of mutations in autosomal inversions")+
  ThemeSobr
MinMutAuto

# Panel e
MeanMutAuto=base+geom_line(data=GoodSimulSubCompleteB[(GoodSimulSubCompleteB$Position=="Autosome"),],
                           aes(x=Gen, y=MeanMutInv, group=Rep), color=Col[1],  size=0.8, alpha=0.5)+
  facet_wrap(.~paste0("U=",U), scales = "free")+
  xlab("Generation")+
  ylab("Number of mutations")+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  labs(title = "Mean number of mutations in autosomal inversions")+
  ThemeSobr
MeanMutAuto

# Panel b
base=ggplot(GoodSimulSubCompleteB[(GoodSimulSubCompleteB$Position=="Y"),])
FreqY=base+geom_line(aes(x=Gen, y=Freq, group=Rep), color=Col[2],  size=0.8, alpha=0.5)+
  facet_wrap(.~paste0("U=",U), scales = "free")+
  xlab("Generation")+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  ylab("Y-linked inversion frequency")+
  labs(title = "Inversion frequency")+
  ThemeSobr
# Panel d
MinMutY=base+geom_line(data=GoodSimulSubCompleteB[(GoodSimulSubCompleteB$Position=="Y" ),],
                       aes(x=Gen, y=MinMutInv, group=Rep),  size=0.8, alpha=0.5, color=Col[2])+
  facet_wrap(.~paste0("U=",U), scales = "free")+
  xlab("Generation")+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  ylab("Number of mutations")+
  labs(title = "Minimum number of mutations in Y-linked inversions")+
  ThemeSobr

# Panel f
MeanMutY=base+geom_line(data=GoodSimulSubCompleteB[(GoodSimulSubCompleteB$Position=="Y"),],
                        aes(x=Gen, y=MeanMutInv, group=Rep), color=Col[2],  size=0.8, alpha=0.5)+
  facet_wrap(.~paste0("U=",U), scales = "free")+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),)+
  xlab("Generation")+
  ylab("Number of mutations")+
  labs(title = "Mean number of mutations in Y-linked inversions")+
  ThemeSobr

# Panel g
u=1e-09
n=2000000
s=0.01
h=0.1
q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
NmutAUT=n*q
t=seq(1,15000,1)
NMut=n*(u/(s*h))*(1-exp(-s*h*t))
DataTableU1=data.frame(t,u,NMut)
u=1e-08
q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
NmutAUT=n*q
t=seq(1,15000,1)
NMut=n*(u/(s*h))*(1-exp(-s*h*t))
DataTableU2=data.frame(t,u,NMut)
DataTable=rbind(DataTableU1,DataTableU2)
NMutDeterministic=ggplot(DataTable)+geom_line(aes(x=t, y=NMut), size=2)+  facet_grid(.~Position)+
  ThemeSobr+
  xlab("Generation")+
  ylab("Number of mutations")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  facet_wrap(.~paste0("u=",u), scales = "free")+labs(title = "Deterministic change in mutation number")

#Merged plot
MergedPlot=plot_grid(FreqAuto, FreqY ,MinMutAuto, MinMutY ,MeanMutAuto, MeanMutY, ncol=2, labels = c('a', 'b', 'c','d','e','f'))
MergedPlot2=plot_grid(MergedPlot, NMutDeterministic, ncol=1, labels=c('','g'),rel_heights = c(2,1),rel_widths = c(2,1))
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS9.png"),MergedPlot2, nrow=4, ncol=2)
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS9.pdf"),MergedPlot2, nrow=4, ncol=2)

#####################################################
################@ Figures S10-13 ####################
#####################################################

Simul=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/X_Y_Autosomal_Inversion_Trajectories_N=1000_FigS10-13.txt",sep=""), stringsAsFactors = F)#Simulation runned with inversion on the X or Y chromosome, or on the autosome. We recorded simulation state every 10 generations.
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv","Chrom", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Simul=Simul[!(Simul$DebInv>10000000 & Simul$Chrom=="Y"),] # 20,000 inversion are runned in autosome (10,000 in Y-bearing genome, 10,000 in X-bearing genomes). Removing here all the simulations that were perform with autosomal inversion on the Y-bearing slim genome. (Same result are obtain by removed autosomal inverison on X-bearing genome)
Simul[Simul$DebInv>10000000,]$Chrom="Autosome" #Label the autosomal inversions
Simul$InvSize=Simul$FinInv - Simul$DebInv #Inversion size
Simul$U=Simul$InvSize*Simul$u
SimulSub=Simul[(Simul$h %in% c(0.5, 0.1, 0.01) & Simul$s %in% c(-0.01, -0.001)),]

SimulSub$Code=paste(SimulSub$N,SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$InvSize,SimulSub$Position, SimulSub$Rep, sep="_") # Define a code identifying eavh simulation
summarySub=SimulSub %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)
FixedSimul=summarySub[summarySub$maxFreq>0.95,]$Code #Grep the code of the inversion that have fixed

LineToAdd2=SimulSub[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate ending states of fixed inversion.

for (i in FixedSimul){ #Do the same thing for fixed inversions
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSub[(SimulSub$Code==i & SimulSub$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimulSub$Gen=24991  #Change it to 24991 (last recorded generation)
  FalseEndGoodSimulSub$Freq=1 #Set its frequency to 1
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub  #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub) 
}

GoodSimulSubComplete=rbind(SimulSub,LineToAdd2) #Add these false simulation end to the simulation data.frame
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(3) #Color Palette
GoodSimulSubComplete[GoodSimulSubComplete$Chromosome=="Autosome",]$Rep=100000+GoodSimulSubComplete[GoodSimulSubComplete$Chromosome=="Autosome",]$Rep

options(scipen=0)
GoodSimulSubComplete$Gen=GoodSimulSubComplete$Gen-15000
DataAllEnd=GoodSimulSubComplete %>% group_by(U, h, s, Rep, Chrom) %>% summarise(MaxGen=max(Gen), MaxFreq=max(Freq))#Grep the end generation of each simulation
DataAllEnd$state="Lost"
DataAllEnd[(DataAllEnd$MaxGen>=9991 & DataAllEnd$MaxFreq==1),]$state="Fixed"
DataAllEnd[(DataAllEnd$MaxGen>=9991 & DataAllEnd$MaxFreq<1),]$state="Segregating"
DataAllEnd[(DataAllEnd$Chrom=="Y" & DataAllEnd$MaxGen>=9991 & DataAllEnd$MaxFreq==0.25),]$state="Fixed"
DataAllEnd[(DataAllEnd$Chrom=="X" & DataAllEnd$MaxGen>=9991 & DataAllEnd$MaxFreq==0.75),]$state="Fixed"

SumEnd=DataAllEnd %>% group_by(U, h, s, Chrom) %>% count(state, sort = TRUE)  #Summarize the data for plotting summary number
SumEnd=SumEnd[SumEnd$state!="Lost",]
SumEnd$Pos=0.10
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chrom=="Y"),]$Pos=0.30 #Define position for plotting the text
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chrom=="Autosome"),]$Pos=1.02
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chrom=="X"),]$Pos=0.77
SumEnd$Parameter=paste0("s=", SumEnd$s, ", h=", SumEnd$h, ", U=", SumEnd$U)
GoodSimulSubComplete$Parameter=paste0("s=", GoodSimulSubComplete$s, ", h=", GoodSimulSubComplete$h, ", U=", GoodSimulSubComplete$U)

### Inversion frequency
for (U in unique(GoodSimulSubComplete$U)){
  
PlotTrajDrift=ggplot(GoodSimulSubComplete[GoodSimulSubComplete$U==U,]) +rasterize(geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chrom)), size=0.2, alpha=0.3), dpi=300)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=1, size=0.2)+
  geom_text(data=SumEnd[SumEnd$U==U,], aes(x=8500, y=Pos, label=paste0("n=",n), color=as.factor(Chrom)), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-linked", "X-linked"))+
  xlab("Generation")+ylab("Inversion frequency")+
  ggtitle("Trajectory of inversions with drift")+
  ThemeSobr+
  theme(
    legend.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title=element_text(size=28, face="bold",hjust=0.5),
    strip.placement = "outside",
    legend.key=element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face="bold", size=20),
    strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
    text = element_text(size=20),
    legend.position="none",
    axis.title = element_text(size = 22, face="bold"))+
  #guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(Parameter ~ Chrom)

save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS10_2_U", U, ".png"), PlotTrajDrift, ncol=3, nrow=6)
}

#####################################################
################@ Figure S14 #########################
#####################################################

Simul10k=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/InversionTrajectories_N=10000_FigS14_2025.txt",sep=""), stringsAsFactors = F)
# Simul10k=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/test.txt",sep=""), stringsAsFactors = F)
colnames(Simul10k)=c("N", "u", "r", "h", "s", "Gen", "StartInv", "EndInv", "Rep", 
                     "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                     "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                     "InvFit", "NoInvFit","Freq","Chromosome")

Simul1k=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/InversionTrajectories_N=1000_Fig3.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
# Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/Compressed/InversionTrajectories_N=1000_Fig3-S13-S15-S19.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000

#Set the column names
colnames(Simul1k)=c("N", "u", "r", "h", "s", "Gen", "StartInv", "EndInv", "Rep", 
                    "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                    "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                    "InvFit", "NoInvFit","Freq","Chromosome")

Simul=full_join(Simul10k, Simul1k)
Simul$Position="Y" #Define position 
Simul[Simul$StartInv>10000000,]$Position="Autosome" 
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

Col5=scales::viridis_pal(begin=0.1, end=0.8, option="A")(5)
Base1k=ggplot(DataSummary[(DataSummary$N==1000),], aes( y=ProbSpread))
# Base1k=ggplot(DataSummary[(DataSummary$s %in% c(-0.001, -0.01, -0.1) & DataSummary$U %in% c(0.001, 0.01, 0.05) & DataSummary$N==1000),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv1k=Base1k+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  geom_point(aes(x=h,color=as.factor(s)), size=5)+
  scale_color_manual("s", values=Col5)+
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

# Base10k=ggplot(DataSummary[(DataSummary$s %in% c(-0.001, -0.01, -0.1) & DataSummary$U %in% c(0.001, 0.01, 0.05) & DataSummary$N==10000),], aes( y=ProbSpread))
Base10k=ggplot(DataSummary[(DataSummary$N==10000),], aes( y=ProbSpread))

PlotPropInvSpread_AllInv10k=Base10k+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  geom_point(aes(x=h,color=as.factor(s)), size=5)+
  scale_color_manual("s", values=Col5)+
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
# 
Merge1k_10k=align_plots(PlotPropInvSpread_AllInv1k, PlotPropInvSpread_AllInv10k + theme(legend.position = "none"), align = 'hv', axis = 'rltp')
PlotMerged=plot_grid(Merge1k_10k[[1]], Merge1k_10k[[2]], legend, rel_widths = c(1,1,0.25), nrow=1)
PlotMerged


title <- ggdraw() +
  draw_label(
    "Fraction of inversions fixed after 10,000 generations",
    fontface = 'bold',
    size=20,
    x = 0,
    hjust = -0.3
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

MergedWtTitle=plot_grid(
  title, PlotMerged,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)
MergedWtTitle

save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS14.pdf"),MergedWtTitle, ncol=2, nrow=3.5)  #Fig. S14
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS14.pdf"),MergedWtTitle, ncol=2, nrow=3.5)  #Fig. S14

#####################################################
################@ Figure S15 #########################
#####################################################

Col=scales::viridis_pal(begin=0, end=0.6, option="A")(2)
Simul=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/RecombinationSuppressorEvolution_n=2Mb_FigS15.txt",sep=""), stringsAsFactors = F) #Evolution of recombination suppressors (suppress recombination in heterozygote and in homozygote, unlike inversions)
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Simul=subset(Simul, Simul$Gen %% 10 == 0) #Keep only every 10 generation for computing purpose
Simul$InvSize=Simul$FinInv - Simul$DebInv #Size of the region affected (n)
SimulSub=Simul #Copy, just in case
SimulSub$Link="Linked" #Define linkage
SimulSub[SimulSub$DebInv > 10000000,]$Link="Unlinked" 
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv, SimulSub$Rep, sep="_") #Define Code
FocS=c(-0.01,-0.05,-0.1) #Value of s to focus on 
SimulSub=SimulSub[(SimulSub$s %in% FocS),] 
SimulSub20G=unique(SimulSub[SimulSub$Gen>15020,]$Code) #Keep only inversion not lost during the first 20 generation
SimulSub=SimulSub[SimulSub$Code %in% SimulSub20G,]
summarySub=SimulSub %>% group_by(Code) %>% summarise(LastFreq=last(Freq), maxGen=max(Gen)) # As before, recreate end state
LostSimulSub=summarySub[summarySub$maxGen<25000,]$Code # Inversion that have been lost or fixed
GoodSimulSub=SimulSub #Keep only non-lost Inversion
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15010),]#For lost inversion, grep their initial state
FalseEndGoodSimulSub$Gen=25000 # Modify their initial state
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
if (length(summarySub[summarySub$LastFreq>0.95,]$Code)>0)
{
  FixedSimul=summarySub[summarySub$LastFreq>0.95,]$Code #Grep the code of the inversion that have fixed
  FalseEndGoodSimulSub[FalseEndGoodSimulSub$Code %in% FixedSimul,]$Freq=1.0 #For all inversion that have fixed, define their end frequency as 1.0
}
GoodSimulSub=rbind(GoodSimulSub,FalseEndGoodSimulSub) #Concatenate the dataset
DataAll=GoodSimulSub
DataAll$s=-DataAll$s
DataAllEnd=subset(DataAll, DataAll$Gen==25000)
DataAllEnd$state="Lost"
DataAllEnd[DataAllEnd$Freq==0.5,]$state="Segregate"
SumEnd=DataAllEnd %>% count(h,s,Link, InvSize, state, sort = TRUE) 
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
  xlab("Generation")+ylab("Inversion frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.1 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.1")
Plot0.1_2Mb

base=ggplot(DataAll[(DataAll$h==0.01 & DataAll$InvSize==2000000),])
Plot0.01_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Inversion frequency")+
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
  xlab("Generation")+ylab("Inversion frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.001 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.001")

Plot2Mb=plot_grid(Plot0.001_2Mb,Plot0.01_2Mb,Plot0.1_2Mb, nrow=2, labels=c("a", "b","c"), label_size = 30)

save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS15.png", Plot2Mb, ncol=4, nrow=6)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS15.pdf", Plot2Mb, ncol=4, nrow=6)




#####################################################
################@ Figure S16 #########################
#####################################################

SimulLamb=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/MutationLambdaDistribution_InversionTrajectories_FigS16.txt",sep=""), stringsAsFactors = F)
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
Base=ggplot(DataSubLamb[DataSubLamb$s<0.5,], aes( y=ProbSpread))
PlotPropInvSpread_AllInvLamb=Base+geom_point(aes(x=s, color=as.factor(U)), size=4)+
  geom_line(aes(x=s, color=as.factor(U)), size=2)+
  scale_color_manual("Region-wide mutation rate (U)", values=Col7)+
  scale_x_log10(breaks=c(0.1, 0.01, 0.001))+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("mean(s)")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid( ~ Position)
PlotPropInvSpread_AllInvLamb
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS16.png"),PlotPropInvSpread_AllInvLamb, ncol=2.3, nrow=3) #Fig. S16
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS16.pdf"),PlotPropInvSpread_AllInvLamb, ncol=2.3, nrow=3) #Fig. S16



#####################################################
################@ Figure S17 #########################
#####################################################

Simul=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/N=1000_2MbChromosomeFusion_Trajectory_FigS17.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
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

save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS17.pdf"),PlotA, nrow=3, ncol=2) #Fig. S17
save_plot(paste0("~/Paper/ModelSexChrom_V2/Plots2/FigS17.png"),PlotA, nrow=3, ncol=2) #Fig. S17


#####################################################
################@ Figure S18-19 #####################
#####################################################

ResultHeatLessLoad=data.frame(h=double(), s=double(), U=double(), q=double(),
                              nq=double(), HC=double(), ExpectedFreqY=double(), ExpectedFreqAutosome=double(),
                              FracFixedY=double(), FracFixedAuto=double(), TotProb=double())
u=2.2/3100000000
for (U in c(0.001, 0.01))
{
  n=as.integer(U/u)
    for (h in seq(0.001,0.25,0.0025))
    {
      for (s in  seq(0.001,0.1,0.001))
      {
        for (HC in c(0, 0.001,0.005,0.01,0.1))
        {
          q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
          nq=n*q
          CumulProbY=0
          CumulProbAuto=0
          FracFixedY=0
          FracFixedAuto=0
          TotProb=0 #Sum of inversion occurring probability: must equal 1 at the end if all inversion are considered (m in seq(0,n))
          for (m in seq(0,nq)) # Only consider less-loaded inversions
          {
            Pm=dbinom(m,n,q) #Proba of occuring inversion with m mutations
            WNI=(exp(m*q*(2*h*s-s)-h*s*(n*q+m)))*(1-HC)
            WNN=exp(n*q*q*(2*h*s-s)-h*s*(n*q*2))
            WII=exp(-s*m)
            if (WNI > WNN)
            {
              FequilY=1 
              FracFixedY=FracFixedY+Pm
              if(WNI>WII)
              {
                s1=WNI-WII
                s2=WNI-WNN
                FequilAuto=s2/(s1+s2)
              }
              else #WII>WNI
              {
                FequilAuto=1
                FracFixedAuto=FracFixedAuto+Pm
              }
            }
            else #WNN>WNI
            { 
              FequilY=0 
              FequilAuto=0
            }
            ProbEquilY=Pm*FequilY
            ProbEquilAuto=Pm*FequilAuto
            CumulProbY=CumulProbY+ProbEquilY
            CumulProbAuto=CumulProbAuto+ProbEquilAuto
            TotProb=TotProb + Pm
          }
          ResultHeatLessLoad[nrow(ResultHeatLessLoad)+1,]=c(h,s,U, q, n*q,HC, CumulProbY,CumulProbAuto,FracFixedY, FracFixedAuto, TotProb)
        }
    }
  }
}


ResultHeatLessLoad$FracFixedY=ResultHeatLessLoad$FracFixedY/ ResultHeatLessLoad$TotProb
ResultHeatLessLoad$FracFixedAuto=ResultHeatLessLoad$FracFixedAuto/ ResultHeatLessLoad$TotProb
LongTableFracFix=ResultHeatLessLoad %>% pivot_longer(cols=c(FracFixedY, FracFixedAuto), names_to="Position", values_to="FracFixed")
LongTableFracFix[LongTableFracFix$Position=="FracFixedY",]$Position="Y-linked"
LongTableFracFix[LongTableFracFix$Position=="FracFixedAuto",]$Position="Autosomal"

nqLimit=data.frame(h=double(),U=double(),HC=double(), s_nq1=double(),
                   s_nq2=double(), s_nq3=double(), s_nq4=double(), s_nq5=double(), s_nq6=double())

for (h in seq(0.001,0.30,0.0025)) # Depend only on h and mu, and not on s
{
  for (U in c(0.01, 0.001))
  {
    for (HC in c(0, 0.001,0.005,0.01))
    {
      s_nq1=(U)/(h*1)
      s_nq2=(U)/(h*2)
      s_nq3=(U)/(h*3)
      s_nq4=(U)/(h*4)
      s_nq5=(U)/(h*5)
      s_nq6=(U)/(h*6)
      nqLimit[nrow(nqLimit)+1,]=c(h,U,HC,s_nq1,s_nq2,s_nq3,s_nq4,s_nq5,s_nq6)
    }
  }
}

colnames(nqLimit)[4]="nq=1" #Change label for plotting purpose
colnames(nqLimit)[5]="nq=2"
colnames(nqLimit)[6]="nq=3"
colnames(nqLimit)[7]="nq=4"
colnames(nqLimit)[8]="nq=5"
colnames(nqLimit)[9]="nq=6"
nqLimitLong=nqLimit %>% pivot_longer(cols=c("nq=1","nq=2","nq=3","nq=4","nq=5","nq=6"), names_to="nq", values_to="Values")

base=ggplot(LongTableFracFix[(LongTableFracFix$HC %in% c(0, 0.001,0.005, 0.01) & LongTableFracFix$U==0.01),])
InversionFrac=base+  geom_tile(aes(x=h, y=s, fill=FracFixed))+
  geom_step(data=nqLimitLong[nqLimitLong$U==0.01,], aes(x=h, y=Values, color=nq), direction='mid', alpha=0.5)+
  scale_color_manual(values=rep("red",6), guide=F)+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  facet_grid(HC~Position)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  scale_fill_viridis( option="D","", direction = -1)+
  labs(title="Fraction of less-loaded inversions expected to fix")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.04),
        legend.key.size = unit(1.2, 'cm'),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold", size=16),
        strip.text = element_text(face="bold", size=18),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        axis.title = element_text(face="bold", size=16),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=18, face="bold",hjust=0.5, vjust=9),
        plot.margin = margin(35, 3, 3, 3, "pt"))+
  guides(fill = guide_colorbar( title = NULL,
                                label.position = "top",
                                label.vjust =1
  ))
InversionFrac
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS18.pdf", InversionFrac, nrow=4, ncol = 2)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS18.png", InversionFrac, nrow=4, ncol = 2)

base=ggplot(LongTableFracFix[(LongTableFracFix$HC %in% c(0, 0.001,0.005, 0.01) & LongTableFracFix$U==0.001),])
InversionFrac=base+ geom_tile(aes(x=h, y=s, fill=FracFixed))+
  geom_step(data=nqLimitLong[nqLimitLong$U==0.001,], aes(x=h, y=Values, color=nq), direction='mid', alpha=0.5)+
  scale_color_manual(values=rep("red",6), guide=F)+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  facet_grid(HC~Position)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  scale_fill_viridis( option="D","", direction = -1)+
  labs(title="Fraction of less-loaded inversions expected to fix")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.04),
        legend.key.size = unit(1.2, 'cm'),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold", size=16),
        strip.text = element_text(face="bold", size=18),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        axis.title = element_text(face="bold", size=16),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=18, face="bold",hjust=0.5, vjust=9),
        plot.margin = margin(35, 3, 3, 3, "pt"))+
  guides(fill = guide_colorbar( title = NULL,
                                label.position = "top",
                                label.vjust =1
  ))
InversionFrac
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS19.pdf", InversionFrac, nrow=4, ncol = 2)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS19.png", InversionFrac, nrow=4, ncol = 2)

#####################################################
################@ Figure S20 ########################
#####################################################

Simul=read.table(paste("~/Paper/ModelSexChrom_V2/Datasets/HaploDiplo_InvTrajectories_FigS20.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000

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


DataSummary=Summary %>% group_by(N,u,r,h,s,InvSize, FreqHaplo, Position) %>% summarise(ProbSpread=mean(StateCode)) # For each set of parameter, compute the fraction of mutation fixed (only for not mutation-free inversion)
options(scipen=0) #Scientific notation
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(4) #Define color

Base=ggplot(DataSummary, aes( y=ProbSpread)) #Plot the result
PlotPropInvSpread_HD=Base+geom_point(aes(x=as.factor(h), color=as.factor(FreqHaplo)), size=3, alpha=0.7)+
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
PlotPropInvSpread_HD

save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS20.png", PlotPropInvSpread_HD, ncol=2, base_aspect_ratio = 1)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS20.pdf", PlotPropInvSpread_HD, ncol=2, base_aspect_ratio = 1)

#####################################################
################@ Figure S21 ########################
#####################################################

data=read.table("~/Paper/ModelSexChrom_V2/Datasets/HaploDiplo_InitStat_FigS21.txt", stringsAsFactors = F)
colnames(data)=c("N", "mu","r","h","s","FreqHaplo", "generation", "meanNbMut", "Nmut1", "MeanFreq1")
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(4)
options(scipen=0)
base=ggplot(data)
nMut=base+geom_line(aes(x=generation, y=Nmut1, color=as.factor(FreqHaplo)))+
  facet_grid(paste0("h=",h)~paste0("s=",s), scale="free")+
  scale_color_manual("Frequency of haploid phases \n (every x generation)", values=Col)+
  ylab("Number of mutations \n segregating genome-wide")+
  xlab("Generation (burn-in)")+
  ThemeSobr

Mean_nMut=base+geom_line(aes(x=generation, y=meanNbMut, color=as.factor(FreqHaplo)))+
  facet_grid(paste0("h=",h)~paste0("s=",s), scale="free")+
  scale_color_manual("Frequency of haploid phases \n (every x generation)", values=Col)+
  ylab("Average number of mutations \n per genome")+
  xlab("Generation (burn-in)")+
  ThemeSobr

FreqMut=base+geom_line(aes(x=generation, y=MeanFreq1, color=as.factor(FreqHaplo)))+
  facet_grid(paste0("h=",h)~paste0("s=",s), scale="free")+
  scale_color_manual("Frequency of haploid phases \n (every x generation)", values=Col)+
  ylab("Mean frequency of mutations")+
  xlab("Generation (burn-in)")+
  ThemeSobr

PlotMerg=plot_grid(nMut, Mean_nMut, FreqMut,  ncol=1, labels = c('a', 'b', 'c'))
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS21.pdf", nrow = 3, ncol=2,PlotMerg)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS21.png", nrow = 3, ncol=2,PlotMerg)

#####################################################
################@ Figure S22 ########################
#####################################################

Dir="~/Paper/ModelSexChrom_V2/Datasets/"
# ID="1750881723915" #Seed Id of the simulation to plot (figure 4) #Used to easily plot multiple simulation
ID="1953648752606" #Seed Id of the simulation to plot (uncomment to create figure S22)

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
save_plot(paste("~/Paper/ModelSexChrom_V2/Plots2/FigS22.png"), MergedPlot, ncol = 2, nrow=2)
save_plot(paste("~/Paper/ModelSexChrom_V2/Plots2/FigS22.pdf"), MergedPlot, ncol = 2, nrow=2)

#####################################################
################@ Figure S23 #########################
#####################################################

data=read.table("~/Paper/ModelSexChrom_V2/Datasets/BurnIn_Stat_FigS23.txt", stringsAsFactors = F) #Stat computed during the burn in of each simulations
colnames(data)=c("N", "mu","r","h","s","generation", "meanNbMut", "Nmut1", "MeanFreq1", "NbMutXY", "FreqMutXY")
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(7)#Color Palette 
options(scipen=0) #Non-scientific notation
base=ggplot(data) #Plot number of mutation
nMut=base+geom_line(aes(x=generation, y=Nmut1, color=as.factor(h)))+ 
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Number of mutations genome-wide")+
  xlab("Generation (burn-in)")+
  ThemeSobr

FreqMut=base+geom_line(aes(x=generation, y=MeanFreq1, color=as.factor(h)))+
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Mean frequency of mutations")+
  xlab("Generation (burn-in)")+
  ThemeSobr

PlotMerg=plot_grid(nMut, FreqMut,  ncol=1, labels = c('a', 'b'))
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS23.png", nrow = 2, ncol=2,PlotMerg)
save_plot("~/Paper/ModelSexChrom_V2/Plots2/FigS23.pdf", nrow = 2, ncol=2,PlotMerg)



