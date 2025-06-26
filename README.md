# MutationShelteringTheoryV2
This is a collection of scripts used to perform the analyses of the manuscript *Stepwise expansion of recombination suppression on sex chromosomes and other supergenes through lower load advantage and deleterious mutation sheltering*.

For any question regarding the analyses or the script usage, please contact me: *paul.yann.jay[at]gmail.com*

For the simulations with SLiM, their is two types of analyses:

# Analyses of the fate of many inversions (Figures 3, S2-3,10-19,22-23,25) 
For these analyses, we first perform burn-in simulations to create "initial populations", and then introduce in the initial populations many inversions or other recombination modifiers.
## Burn In
These initial populations can be created with the scripts:
 - ScriptNeutralInversion_DefineInitialState_XY_ModForSlim4.3.slim #Basic simulation (Figures 2,3, S2-3, S10-13)
 - ScriptNeutralInversion_DefineInitialState_HaploDiplo_SingleChrom.slim #Simulation with HaploDiplontic life cycle 
 - ScriptNeutralInversion_DefineInitialState_TwoChromFus.slim #Simulation of chromosome fusion (Figure S18)
 - ScriptNeutralInversion_DefineInitialState_XY_Variable_sh.slim #Simulation with mutation of variable h and s (drawn from distribution)
 - ScriptNeutralInversion_DefineInitialState.slim #Simulation of populations with non-XY system but with two mating types (all individuals are always heterozygous) (Figure S14)

Basic usage of the "DefineInitialState" script:
>slim -d N=1000 -d s=-0.01 -d h=0.1 -d r=1e-08 -d mu=1.45e-8 ScriptNeutralInversion_DefineInitialState_XY.slim #Perform a single simulation  
parallel -j10 slim -d N=1000 -d mu=1e-8 -d h={1} -d s={2} -d r=1e-6 ScriptNeutralInversion_DefineInitialState_XY.slim ::: 0 0.001 0.005 0.01 0.05 0.1 0.5 ::: -0.001 -0.005 -0.01 -0.05 -0.1 -0.6 #Perform multiple simulation in parallel on 10 cores with different *s* and *h*  

*mu* define the mutation rate, *r* the recombination rate, *s* the selection coefficient (fixed or the mean of the gamma distribution), *h* the dominance coefficient, *N* the population size

This writes the initial populations in a ../InitialState directory (must be created before!)

## Introduce inversions
Then, we introduce in these initial populations 1000's of inversions or other recombination modifiers. To do so, use the scripts:
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_ModForSlim4.3.slim #Basic simulation (Figures 3, S2-3, S10-13 )
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_SingleChrom_HaploDiplo.slim #Simulation with HaploDiplontic life cycle (Figures S22-23)
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_Variablesh.slim #Simulation with mutation of variable h and s (drawn from distribution, Figure S17)
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_ChromFus.slim #Simulation of chromosome fusion (Figure S18)
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_HomoHeteroRecombMod.slim #Simulation of populations with non-XY system but with two mating types (all individuals are always heterozygous) and the introduction of not inversion but recombination modifier suppresssing recombination whatever their genotype (Figure S14)

Basic usage of the "IntroduceInversionFromInitStat" scripts :
 > slim -d N=1000 -d mu=1e-08 -d h=0.01 -d s=-0.01 -d r=1e-5 -d rep=1 -d start=1 -d end=2000000 -d "Init='../InitialState/slim_g15000_N=1000_r=1e-05_u=1e-08_s=-0.01_h=0.01_1635828730194.txt'" -d "SexChrom='Y'" Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY.slim #Introduction of a single 2000000 bp inversion on the Y chromosome in a initial population. Rep is only used to annotate the simulation (simunation number "Rep"). 

The position of the inversion is defined with *start* and *end*. To mimick chromosome fusion, just consider start and end spanning the gap between the two chromosome, *i.e.* the position 10Mb (*e.g. start=7000000, end=11000000)

To introduce 1000's of inversion in population with different parameters, use (and adapt) the script ParalleleWhileLoop_IndivSimulPlot_NMut_XY.sh:
 > source ParalleleWhileLoop_IndivSimulPlot_NMut_XY.sh  
  export -f slimFc  
  parallel -j30 slimFc ::: 1e-08 ::: 0.5 ::: -0.005 -0.001 -0.01 -0.05 -0.1 ::: {1..10000} ::: 5000000 15000000 ::: 500000 2000000 1000000 5000000 ::: Y X ## Run on 30 cores, 10000 inversion per parameter combinations are simulated, with one h and five possible s, inversions of 4 different sizes are considered, either on a X-bearing genome or in a Y-bearing genome, and either in the sex-chromosome (Mid position of the inversion at 5000000) or in the autosome (Mid position of the inversion at 15000000) #It automatically detect the right initial population to use, if it exists.  

This produces output in ../Output (must be created) such as N=1000_Inv=4000001-6000001_r=1e-06_u=1e-09_h=0_s=-0.25IntroduceInvFromInit_Nmut_Freq_IndivPlot.txt #file containing the trajectory of 10,000 2Mb sex-linked inversions with N=1000, s=-0.25, h=0, etc # Result of the Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY.slim script. Other script produce output with slightly different names.
For plotting, concatenate the output
 > cat N=1000_*_IndivPlot.txt > InversionTrajectories_N=1000.txt #Can be plotted with the R script to produce figure 3.
	
# Analyses of the formation of sex chromosomes (Figures 4-5 and S24):
In this case, just run the script ScriptFormationXYChromosome_VarGamma_OnlyXY_Revers_Optimized.slim with parameters. For instance:

 > slim -d N=10000 -d mu=1e-09 -d s=-0.03 -d r=1e-06 -d rep=5 -d N_BP=10000 -d InvR=0.00000001 -d MaxSizeInv=20000000 ScriptFormationXYChromosome_VarGamma_OnlyXY_Revers_Optimized.slim # N_BP is the number of breakpoint in the genome, InvR is the inversion rate, MaxSizeInv is the maximum inversion size. 

The outputs of these simulation need to be parsed for plotting. For that, use the Parse* scripts:
 - ParseInvFreqOutput_Rev.pl
 - ParseRecombinationOutput_SepSex.pl 

such that:
 > ./ParseInvFreqOutput.pl -i N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_InvFreq_IndivSimulation_XY.txt -o N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_InvFreq_IndivSimulation_XY.Parsed.txt  
./ParseRecombinationOutput.pl -i N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.txt -o N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.Parsed.txt

# Deterministic simulations (Figure 2, S1,4-9) and Figures
Figures and deterministic simulations can be performed with the R script (more details in the script):
 - AllFigureV5.R # For the main figures
 - SupFig_V2.R # For the sup figures

The original datasets from Jay et al are available on figshare (doi:10.6084/m9.figshare.19961033). The new datasets will be published upon publication, but by then are available upon request.
