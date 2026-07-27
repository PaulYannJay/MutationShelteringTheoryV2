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

