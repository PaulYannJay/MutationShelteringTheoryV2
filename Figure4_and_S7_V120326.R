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

df$NeutralExpect=1/(2*df$N)
df[df$Chrom=="Y chromosome",]$NeutralExpect=2/df$N

df$NormProbFix=df$PropFix/df$NeutralExpect
#Conf interval
ntrials=df$Fixed+df$NotFixed
df$PropFix_lo=qbeta(0.025, df$Fixed+0.5, df$NotFixed+0.5)
df$PropFix_hi=qbeta(0.975, df$Fixed+0.5, df$NotFixed+0.5)

## Transform to the normalized scale
df$NormProbFix_lo=df$PropFix_lo/df$NeutralExpect
df$NormProbFix_hi=df$PropFix_hi/df$NeutralExpect

df_plot = df[(df$sInv %in% c(0.0, 0.01)), ]

df_plot$s_jit = df_plot$s

df_plot[(df_plot$hLegend == "With fully recessive mutations" & df_plot$s==-0.01 & df_plot$U==0.024), ]$s_jit = -0.0085
df_plot[(df_plot$hLegend == "With fully recessive mutations" & df_plot$s==-0.01 & df_plot$U==0.06), ]$s_jit = -0.0095
df_plot[(df_plot$hLegend == "Without fully recessive mutations" & df_plot$s==-0.01 & df_plot$U==0.024), ]$s_jit = -0.0105
df_plot[(df_plot$hLegend == "Without fully recessive mutations" & df_plot$s==-0.01 & df_plot$U==0.06), ]$s_jit = -0.0115
df_plot[(df_plot$hLegend == "With fully recessive mutations" & df_plot$s==-0.001 & df_plot$U==0.024), ]$s_jit = -0.00025
df_plot[(df_plot$hLegend == "With fully recessive mutations" & df_plot$s==-0.001 & df_plot$U==0.06), ]$s_jit = -0.00075
df_plot[(df_plot$hLegend == "Without fully recessive mutations" & df_plot$s==-0.001 & df_plot$U==0.024), ]$s_jit = -0.00135
df_plot[(df_plot$hLegend == "Without fully recessive mutations" & df_plot$s==-0.001 & df_plot$U==0.06), ]$s_jit = -0.00185

library(colorspace)
Col2=c("#f4a261","#2a9d8f")
Col2_light = lighten(Col2, amount = 0.4)
base = ggplot(
  df_plot,
  aes(
    x = s_jit,
    y = NormProbFix,
    fill = as.factor(hLegend),
    shape = as.factor(U)
  )
)
v

PlotProbFix2 = base +
  #geom_hline(aes(yintercept = 1, linetype="Neutral expectation"), size=0.3)+
	geom_errorbar(
		aes(ymin = NormProbFix_lo, ymax = NormProbFix_hi,
		    color=hLegend),
		width = 0,
		linewidth = 0.7
	) +
	geom_point(size = 8, colour = "grey10", alpha=0.9) +
	scale_fill_manual("Distribution of dominance coefficients", values = Col2) +
	scale_color_manual("Distribution of dominance coefficients", values = Col2_light) +
	scale_shape_manual("U", values = c(21, 24)) +
  scale_linetype_manual("",values="dashed")+
	expand_limits(y = 0) +
  scale_x_continuous(
    breaks = c(-0.001, -0.01),
    labels = c("-0.001", "-0.01"),
    expand = expansion(mult = c(0.3, 0.3))
  )+
	ylab("Normalised inversion fixation probability") +
	xlab("Selection coefficient of mutations (s)") +
	themeSimple +
	theme(
		panel.spacing = unit(1.5, "lines"),
		legend.title = element_text(face = "bold"),
		strip.text = element_text(face = "bold"),
		strip.background = element_rect(linewidth = 1, fill = "transparent", color = "black")
	) +
	facet_grid(sInvClean ~ Status, scales = "free")

PlotProbFix2


save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/Figure3Sup.pdf"), PlotProbFix2, nrow=2, ncol=2)    


base = ggplot(
  df_plot[df_plot$N==12500,],
  aes(
    x = s_jit,
    y = NormProbFix,
    fill = as.factor(hLegend),
    shape = as.factor(U)
  )
)


PlotProbFixMain = base +
  #geom_hline(aes(yintercept = 1, linetype="Neutral expectation"), size=0.3)+
  geom_errorbar(
    aes(ymin = NormProbFix_lo, ymax = NormProbFix_hi,
        color=hLegend),
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(size = 8, colour = "grey10", alpha=0.9) +
  scale_fill_manual("Distribution of dominance coefficients", values = Col2) +
  scale_color_manual("Distribution of dominance coefficients", values = Col2_light) +
  scale_shape_manual("U", values = c(21, 24)) +
  scale_linetype_manual("",values="dashed")+
  expand_limits(y = 0) +
  scale_x_continuous(
    breaks = c(-0.001, -0.01),
    labels = c("-0.001", "-0.01"),
    expand = expansion(mult = c(0.3, 0.3))
  )+
  ylab("Normalised inversion fixation probability") +
  xlab("Selection coefficient of mutations (s)") +
  themeSimple +
  theme(
    panel.spacing = unit(1.5, "lines"),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(linewidth = 1, fill = "transparent", color = "black")
  ) +
  facet_grid(sInvClean ~ Chrom, scales = "free")


save_plot(paste0("~/Project/MutationShelteringV2/Revision1/Plots/Figure4.pdf"), PlotProbFixMain, nrow=2, ncol=2)    
