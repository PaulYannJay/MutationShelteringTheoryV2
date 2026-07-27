
### V07/05/26 
library(cowplot)
library(ggplot2)
library(tidyverse)
library(viridis)
library(directlabels)
library(Hmisc)
library(patchwork)
library(scales)
library(ggh4x)
library(grid)

options(scipen = 0)
ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12),
  axis.line = element_line(colour = "grey"),
  legend.spacing.y= unit(0, 'cm'),
  strip.background = element_rect(linewidth = 1,fill='transparent', color="black"),
  strip.text.x = element_text(face="bold", size=16),
  strip.text.y = element_text( size=16),
  #panel.spacing = unit(0.8, "lines"),
  axis.title.y = element_text(face="bold", size=14),
  axis.title.x = element_text(size=14),
  legend.title = element_text(size = 14, face="bold"),
  legend.text = element_text(size = 14)
)
Col3=scales::viridis_pal(begin=0.1, end=0.8, option="A")(3)

#Simul=read.table(paste("~/Project/MutationShelteringV2/Revision1/Output/IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat.24Col.txt",sep=""), stringsAsFactors = F) # Data for N=10000
Simul=read.table(paste("~/Project/MutationShelteringV2/Revision1/Output/files_N10000_s-0.0001/IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat_Full.24Col.txt",sep=""), stringsAsFactors = F) # Data for N=10000
colnames(Simul)=c("N", "u", "r", "h", "s", "sInv", "Gen", "StartInv", "EndInv", "Rep" ,"Chromosome",
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")
SimulCp=Simul # Just in case something wrong happen, to not load again the file
Simul$Gen=Simul$Gen - 6*Simul$N
Simul[Simul$Chromosome=="SexChrom",]$Freq=Simul[Simul$Chromosome=="SexChrom",]$Freq * 4 # Define Y inversion frequency
Simul$InvSize=Simul$EndInv - Simul$StartInv #Inversion size
Simul[Simul$Chromosome=="SexChrom",]$Chromosome="Ychromosome"

#To create FigureS4, uncomment this line:
#Simul=Simul[!(Simul$Chromosome=="Ychromosome" & Simul$Gen > (Simul$N * 6)/4),] #Remove generation after 6N/4 for Y chromosome; Remove this line to consider fixation time probability after 6N generation on both Y chromosome and autosome


Summary=Simul %>% group_by(N,u,r,h,s,InvSize,Chromosome, Rep) %>% summarise(lastFreq=last(Freq), maxGen=last(Gen))
Summary$State="Lost"
Summary[Summary$maxGen==6*Summary$N,]$State="Segregating"
Summary[(Summary$Chromosome=="Ychromosome" & Summary$maxGen==6*Summary$N/4),]$State="Segregating" #Remove this line to consider fixation time probability after 6N generation on both Y chromosome and autosome
Summary[Summary$lastFreq>0.95,]$State="Fixed"
Summary$StateCode=1
Summary[Summary$State=="Lost",]$StateCode=0
Summary[Summary$State=="Segregating",]$StateCode=0
Summary$Norms=Summary$N*Summary$s
Summary=Summary[(Summary$Norms %in% c(-1, -10,-100) & Summary$InvSize %in% c(1000000, 5000000)),]

SegregatingAuto=Summary[(Summary$State=="Segregating" & Summary$Chromosome=="Autosome"),]
SegregatingY=Summary[(Summary$State=="Segregating" & Summary$Chromosome=="Ychromosome"),]

SegregatingLowh=Summary[(Summary$State=="Segregating" & Summary$N==10000),]
Fixed=Summary[Summary$State=="Fixed",]
FixedY=Summary[(Summary$State=="Fixed" & Summary$Chromosome=="Ychromosome"),]
FixedAuto=Summary[(Summary$State=="Fixed" & Summary$Chromosome=="Autosome"),]

DataSummary= Summary %>% group_by(N,u,r,h,s,InvSize, Chromosome) %>% summarise(ProbSpread=mean(StateCode), n=n(),InvFixed=sum(StateCode))
DataSummary$Norms=DataSummary$N*DataSummary$s
DataSummary=DataSummary[DataSummary$Norms %in% c(-1, -10,-100),]
mean(DataSummary[(DataSummary$Chromosome=="Ychromosome"),]$ProbSpread)/mean(DataSummary[(DataSummary$Chromosome=="Autosome"),]$ProbSpread) 
mean(DataSummary[(DataSummary$Chromosome=="Ychromosome" & DataSummary$h<=0.1),]$ProbSpread)/mean(DataSummary[(DataSummary$Chromosome=="Autosome"& DataSummary$h<=0.1),]$ProbSpread) 

mean(DataSummary[(DataSummary$N==1000 & DataSummary$Chromosome=="Ychromosome"),]$ProbSpread)/mean(DataSummary[(DataSummary$N==1000 & DataSummary$Chromosome=="Autosome"),]$ProbSpread) #On average, inversions were 18 time more likely to fix on Y chromosomes than on autosomes with N=1000 and 128 times more likely when N=10,000 
mean(DataSummary[(DataSummary$N==10000 & DataSummary$Chromosome=="Ychromosome"),]$ProbSpread)/mean(DataSummary[(DataSummary$N==10000 & DataSummary$Chromosome=="Autosome"),]$ProbSpread) # On average, inversions were 18 time more likely to fix on Y chromosomes than on autosomes with N=1000 and 128 times more likely when N=10,000 

mean(DataSummary[(DataSummary$N==1000 & DataSummary$Chromosome=="Ychromosome" & DataSummary$h<=0.1),]$ProbSpread)/mean(DataSummary[(DataSummary$N==1000 & DataSummary$Chromosome=="Autosome" & DataSummary$h<=0.1),]$ProbSpread) #On average, inversions were 18 time more likely to fix on Y chromosomes than on autosomes with N=1000 and 128 times more likely when N=10,000 
mean(DataSummary[(DataSummary$N==10000 & DataSummary$Chromosome=="Ychromosome" & DataSummary$h<=0.1),]$ProbSpread)/mean(DataSummary[(DataSummary$N==10000 & DataSummary$Chromosome=="Autosome" & DataSummary$h<=0.1),]$ProbSpread) # On average, inversions were 18 time more likely to fix on Y chromosomes than on autosomes with N=1000 and 128 times more likely when N=10,000 

DataSummary$LowInterval=binconf(DataSummary$ProbSpread*DataSummary$n, n=DataSummary$n,
                                alpha=.05, include.x = F)[,2] #Binomial Confidence interval
DataSummary$HighInterval=binconf(DataSummary$ProbSpread*DataSummary$n, n=DataSummary$n,
                                 alpha=.05, include.x = F)[,3] #Binomial Confidence interval
DataSummary$NeutralFixProb=1/(2*DataSummary$N)
DataSummary[DataSummary$Chromosome=="Ychromosome",]$NeutralFixProb=1/(DataSummary[DataSummary$Chromosome=="Ychromosome",]$N/2)
DataSummary$RelatNeutralFix=DataSummary$ProbSpread/DataSummary$NeutralFixProb #Get the normalised fixation probability by dividing fixation probabilities by neutral expectation
DataSummary$RelatNeutralFixLow=DataSummary$LowInterval/DataSummary$NeutralFixProb
DataSummary$RelatNeutralFixHigh=DataSummary$HighInterval/DataSummary$NeutralFixProb
DataSummary$U=DataSummary$u * DataSummary$InvSize #Region-Wide mutation rate
DataSummary[DataSummary$ProbSpread==0,]$RelatNeutralFixHigh=0
DataSummary[DataSummary$ProbSpread==0,]$HighInterval=0
DataSummary$UL=paste0("U[L]==",DataSummary$U)


mysqrt_trans <- function()  {
  trans_new(
    "mysqrt",
    transform = base::sqrt,
    inverse = function(x) ifelse(x < 0, 0, x^2),
    domain = c(0, Inf)
  )
}


DataSummary$Norms=factor(DataSummary$Norms, levels=c(-1, -10, -100))

Base1k=ggplot(DataSummary[(DataSummary$Norms %in% c(-100, -10, -1) & DataSummary$U %in% c(0.001, 0.01, 0.05) & DataSummary$N==1000),], aes( y=ProbSpread))
PlotProbInvSpread_AllInv1k=Base1k+
  geom_line(data=DataSummary[DataSummary$N==1000,], aes(x=h, y=NeutralFixProb), linetype="dashed")+
  geom_hline(yintercept=0, size=0.4, color="grey")+
  geom_line(aes(x=h, color=as.factor(Norms)), size=0.8)+
  geom_point(aes(x=h,color=as.factor(Norms)), size=4)+
  geom_errorbar(aes(x=h, ymin=LowInterval, ymax=HighInterval, color=as.factor(Norms)), width = 0)+
  scale_color_manual("Ns", values=Col3)+
  ylab("Fraction of inversions fixed after 10,000 generations (square-root scale)")+
  xlab(expression(paste("Dominance coefficient (",bolditalic("h"), ")")))+
  ThemeSobr+
  theme(legend.position = "none", 
        strip.text.y = element_blank(),
        axis.title.x = element_text(face="bold"), 
        panel.grid.major.y = element_line(size=0.2, color="grey90"),
        #panel.background = element_rect(fill="grey92"),
        plot.title = element_text(face = "bold", size=18), 
  )+
  facet_grid(UL ~ Chromosome, scale="free", labeller = labeller(UL=label_parsed, Chromosome = c(Ychromosome = "Y chromosome")))+
  ggtitle("N=1000") +  facetted_pos_scales(
    y = list(
      scale_y_continuous(trans = mysqrt_trans(), limits = c(0, NA), breaks = c(0.0, .0003, .001, .002, .003)),
      scale_y_continuous(trans = mysqrt_trans(), limits = c(0, NA), breaks = c(0,0.0015, 0.005, 0.01, 0.015)),
      scale_y_continuous(trans = mysqrt_trans(), limits = c(0, NA), breaks = c(0, 0.003, 0.01, 0.02, 0.03))
    )
  )
#PlotProbInvSpread_AllInv1k

Base10k=ggplot(DataSummary[(DataSummary$Norms %in% c(-100, -10, -1) & DataSummary$U %in% c(0.001, 0.01, 0.05) & DataSummary$N==10000),], aes( y=ProbSpread))
PlotProbInvSpread_AllInv10k=Base10k+
  geom_line(data=DataSummary[DataSummary$N==10000,], aes(x=h, y=NeutralFixProb, linetype = "Neutral expectation"))+
  geom_hline(yintercept=0, size=0.4, color="grey")+
  geom_line(aes(x=h, color=as.factor(Norms)), size=0.8)+
  geom_point(aes(x=h,color=as.factor(Norms)), size=4)+
  geom_errorbar(aes(x=h, ymin=LowInterval, ymax=HighInterval, color=as.factor(Norms)), width = 0)+
  #geom_hline(yintercept = 1, linetype="dashed", size=0.1)+
  scale_color_manual("Ns", values=Col3)+
  scale_linetype_manual(name = NULL,values = c("Neutral expectation" = "dashed"))+
  ylab("Fraction of inversions fixed after 10,000 generations (square-root scale)")+
  xlab(expression(paste("Dominance coefficient (",bolditalic("h"), ")")))+
  ThemeSobr+
  #scale_y_sqrt(expand(c(0,0)))+
  theme(legend.position = "none", 
        #axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold"), 
        panel.grid.major.y = element_line(size=0.2, color="grey90"),
        plot.title = element_text(face = "bold", size=18))+
  facet_grid(UL ~ Chromosome, scale="free", labeller = labeller(UL=label_parsed, Chromosome = c(Ychromosome = "Y chromosome")))+
  #facet_grid(UL ~ Chromosome, scale="free", labeller=label_parsed)+
  ggtitle("N=10,000") + facetted_pos_scales(
    y = list(
      scale_y_continuous(trans = mysqrt_trans(), limits = c(0, NA), breaks = c(0.0, .0003, .001, .002)),
      scale_y_continuous(trans = mysqrt_trans(), limits = c(0, NA), breaks = c(0,0.0006,0.002, 0.004, 0.006)),
      scale_y_continuous(trans = mysqrt_trans(), limits = c(0, NA), breaks = c(0, 0.0015,0.005, 0.01, 0.015))
      # etc, one per row/col facet depending on your layout
    )
  )

#PlotProbInvSpread_AllInv10k


DataSummaryWide2= DataSummary%>% pivot_wider(id_cols = c(N,u,r,h,s,U,UL,n), names_from = Chromosome, values_from = InvFixed)

mc_log2fc_norm <- function(kY, nY, kA, nA, qY, qA, nsim = 50000) { #function used to calculate the confidence interval of Y/Autosome fixation probability ratio
  
  if (kY == 0 && kA == 0) {
    return(rep(0, nsim))
  }
  
  pY <- rbeta(nsim, kY + 0.5, nY - kY + 0.5)
  pA <- rbeta(nsim, kA + 0.5, nA - kA + 0.5)
  
  log2((pY / qY) / (pA / qA))
}

DataSummaryWide2=DataSummaryWide2[!is.na(DataSummaryWide2$Ychromosome),]
DataSummaryWide2 <- DataSummaryWide2 %>%
  rowwise() %>%
  mutate(
    sims = list(mc_log2fc_norm(Ychromosome, n, Autosome, n, 1/(N/2)*n, 1/(2*N)*n)),
    log2FC = median(sims, na.rm=T),
    lo = quantile(sims, 0.025, na.rm=T),
    hi = quantile(sims, 0.975, na.rm=T)
  ) %>%
  ungroup() %>% select(-sims)

set.seed(10)  # for reproducibility

DataSummaryWide2$Norms=DataSummaryWide2$s*DataSummaryWide2$N
DataSummaryWide2$Norms=factor(DataSummaryWide2$Norms, levels=c(-1, -10, -100, -1000))

DataSummaryWide2$h_jit <- DataSummaryWide2$h 
DataSummaryWide2[DataSummaryWide2$Norms==-1,]$h_jit <- DataSummaryWide2[DataSummaryWide2$Norms==-1,]$h + 0.01
DataSummaryWide2[DataSummaryWide2$Norms==-10,]$h_jit <- DataSummaryWide2[DataSummaryWide2$Norms==-10,]$h - 0.01

signed_sqrt_trans <- function() { #Transformation to apply to scale
  trans_new(
    "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse   = function(x) sign(x) * (x^2),
    domain = c(-Inf, Inf)
  )
}

panel_limits1k <- DataSummaryWide2 %>% #Define the panel limits for y axis
  filter(N == 1000, Norms %in% c(-100,-10,-1)) %>%
  group_by(UL) %>%
  summarise(
    ymin = min(lo, na.rm = TRUE),
    ymax = max(hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(M = ceiling(1.1 * pmax(abs(ymin), abs(ymax))))
  #mutate(M = ceiling(pmax(abs(ymin), abs(ymax))))

y_scales1k <- lapply(seq_len(nrow(panel_limits1k)), function(i) {  #Define the scale for y axis
  ULval <- panel_limits1k$UL[i]
  Mval  <- panel_limits1k$M[i]
  
  as.formula(sprintf(
    "UL == '%s' ~ scale_y_continuous(trans = signed_sqrt_trans(), limits = c(-%s, %s))",
    ULval, Mval, Mval
  ))
})

panel_limits10k <- DataSummaryWide2 %>% #Define the panel limits for y axis
  filter(N == 10000, Norms %in% c(-100,-10,-1)) %>%
  group_by(UL) %>%
  summarise(
    ymin = min(lo, na.rm = TRUE),
    ymax = max(hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(M = ceiling(1.1 * pmax(abs(ymin), abs(ymax))))
  #mutate(M = ceiling(pmax(abs(ymin), abs(ymax))))

y_scales10k <- lapply(seq_len(nrow(panel_limits10k)), function(i) { #Define the scale for y axis
  ULval <- panel_limits10k$UL[i]
  Mval  <- panel_limits10k$M[i]
  
  as.formula(sprintf(
    "UL == '%s' ~ scale_y_continuous(trans = signed_sqrt_trans(), limits = c(-%s, %s))",
    ULval, Mval, Mval
  ))
})

DataSummaryWide2 <- DataSummaryWide2 %>%
  mutate(
    zero_case = case_when(
      Ychromosome > 0 & Autosome == 0 ~ "Positive",
      Ychromosome == 0 & Autosome > 0 ~ "Negative",
      TRUE ~ NA_character_
    )
  )

DataSummaryWide2 <- DataSummaryWide2 %>% #Make the diamond on plots when one of the two region has null fixation probabilities
  mutate(
    star_y = case_when(
      zero_case == "Positive" ~ hi + 0.2 * pmax(abs(lo), abs(hi), 1),
      zero_case == "Negative" ~ lo - 0.2 * pmax(abs(lo), abs(hi), 1))
  )


Col4=scales::viridis_pal(begin=0.1, end=0.8, option="A")(4)
Base1k <- ggplot(DataSummaryWide2[(DataSummaryWide2$Norms %in% c(-100,-10,-1) & DataSummaryWide2$N==1000),], aes(y = log2FC, color = Norms))
PlotPropInvSpreadNeutralCompar_AllInv1k <- Base1k +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2, fill = "#c1121f", alpha = 0.1) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#003049", alpha = 0.1) +
  geom_hline(aes(yintercept = 0, linetype = "Neutral expectation")) +
  scale_linetype_manual(name = NULL, values = c("Neutral expectation" = "dashed")) +
  geom_line(aes(x = h_jit, group = Norms), size = 0.8) +
  geom_point(aes(x = h_jit), size = 4) +
  geom_errorbar(aes(x = h_jit, ymin = lo, ymax = hi), width = 0, size = 0.6, alpha = 0.5) +
  geom_point(
    data = DataSummaryWide2 %>%
      filter(Norms %in% c(-100, -10, -1), N == 1000, !is.na(zero_case)),
    aes(x = h_jit, y = star_y, fill = zero_case, color=Norms),
    shape = 9,
    size = 1.5,
    stroke = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_manual("Ns", values = Col3) +
  scale_fill_manual(
    "Zero count used",
    values = c("Y = 0" = "goldenrod2", "A = 0" = "firebrick2", "both zero" = "black")
  ) +
  ylab(expression(bold(paste(Log[2], " fold change (Y/Autosome)")))) +
  xlab(expression(paste("Dominance coefficient (", bolditalic("h"), ")"))) +
  ThemeSobr +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(size = 0.3, color = "grey95")
  ) +
  facet_grid2(
    UL ~ "Normalised comparisons",
    scales = "free_y",
    independent = "y",
    labeller = labeller(UL = label_parsed)
  ) +
  facetted_pos_scales(y = y_scales1k)


Base10k <- ggplot(DataSummaryWide2[(DataSummaryWide2$Norms %in% c(-100,-10,-1) & DataSummaryWide2$N==10000),], aes(y = log2FC, color = Norms))
PlotPropInvSpreadNeutralCompar_AllInv10k <- Base10k +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-2, fill="#c1121f", alpha=0.1) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,   ymax=Inf, fill="#003049", alpha=0.1) +
  geom_hline(aes(yintercept=0, linetype="Neutral expectation"))+
  scale_linetype_manual(name = NULL,values = c("Neutral expectation" = "dashed"))+
  geom_line(aes(x = h_jit, group = Norms), size = 0.8) +
  geom_point(aes(x = h_jit), size = 4) +
  geom_errorbar(aes(x = h_jit, ymin = lo, ymax = hi), width = 0, size = 0.6, alpha=0.5) +
  geom_point(
    data = DataSummaryWide2 %>%
      filter(Norms %in% c(-100, -10, -1), N == 10000, !is.na(zero_case)),
    aes(x = h_jit, y = star_y, fill = zero_case, color=Norms),
    shape = 9,
    size = 1.5,
    stroke = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_manual("Ns", values = Col3) +
  scale_fill_manual("Ns", values = Col3) +
  ylab(expression(bold(paste(Log[2], " fold change (Y/Autosome)")))) +
  xlab(expression(paste("Dominance coefficient (", bolditalic("h"), ")"))) +
  ThemeSobr +
  theme(        
                legend.position = "none",
                axis.title.x = element_text(face="bold"), 
                axis.title.y=element_blank(),
                strip.text.x = element_blank(),
                panel.grid.major.y = element_line(size=0.3, color="grey95")) +
  facet_grid2(
    UL ~  "Normalised comparisons",
    scales = "free_y",
    independent = "y",
    labeller = labeller(
      UL = label_parsed,
    )) +
  facetted_pos_scales(y = y_scales10k)
PlotPropInvSpreadNeutralCompar_AllInv10k


PlotProbInvSpread_AllInv1k  = PlotProbInvSpread_AllInv1k  + theme(axis.title.y = element_blank())
PlotProbInvSpread_AllInv10k = PlotProbInvSpread_AllInv10k + theme(axis.title.y = element_blank())
PlotPropInvSpreadNeutralCompar_AllInv1k  = PlotPropInvSpreadNeutralCompar_AllInv1k + theme(axis.title.y = element_blank())
PlotPropInvSpreadNeutralCompar_AllInv10k = PlotPropInvSpreadNeutralCompar_AllInv10k  + theme(axis.title.y = element_blank())

# Make a vertical y-title grob
ytitle_leftcol = wrap_elements(
  grid::textGrob(
    "Fraction of inversions fixed after 10,000 generations (square-root scale)",
    rot = 90,
    gp = gpar(fontsize = 16, fontface="bold")
  )
)

ytitle_midcol = wrap_elements(
  grid::textGrob(
    expression(bold(Log[2]*" fold change (Y/Autosome, square-root scale)")),
    rot = 90,
    gp = gpar(fontsize = 16)
  )
)


layout2=c(
  area(1,1,2,1),
  area(1,2,1,2),
  area(1,3,2,3),
  area(1,4,1,4),
  area(2,2,2,2),
  area(2,4,2,4)
)

PlotMerged3=ytitle_leftcol+
  PlotProbInvSpread_AllInv1k+
  ytitle_midcol+
  PlotPropInvSpreadNeutralCompar_AllInv1k+
  PlotProbInvSpread_AllInv10k+
  PlotPropInvSpreadNeutralCompar_AllInv10k+
  plot_layout(design = layout2, widths = c(0.08, 1, 0.08, 0.7)) + 
  plot_annotation(tag_levels = list(c("", "a", "", "b", "","")))&
  theme(plot.tag = element_text(size = 16, face="bold"))
PlotMerged3

save_plot("~/Project/MutationShelteringV2/Revision1/Plots/Merged3Figure2_050526.png", PlotMerged3, ncol=0.8*3, nrow=0.8*4)

save_plot("~/Project/MutationShelteringV2/Revision1/Plots/Merged3Figure2_050526.pdf", PlotMerged3, ncol=0.8*3, nrow=0.8*4)




### With infinite stars ###
