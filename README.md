# MutationShelteringTheoryV2
This is a collection of scripts used to perform the analyses of the manuscript *Stepwise expansion of recombination suppression on sex chromosomes and other supergenes through lower load advantage and deleterious mutation sheltering*.

For any question regarding the analyses or the script usage, please contact me: *paul.yann.jay[at]gmail.com*

For the simulations with SLiM, their is two types of analyses:

# Analyses of the fate of many inversions (Figures 3, S2-3,10-19,22-23,25) 
For these analyses, we first perform burn-in simulations to create "initial populations" under different parameter values, and then introduce in the initial populations many inversions or other recombination modifiers.

## Burn In
These initial populations can be created with the scripts:
 - ScriptNeutralInversion_DefineInitialAutosome_AllDistrib.slim #Basic simulation 

 Basic usage of the "DefineInitialState" script:
> slim -d N=1000 -d s=-0.01 -d h=0.1 -d r=1e-08 -d mu=1.45e-8 ScriptNeutralInversion_DefineInitialStateAutosome_AllDistrib.slim  #Perform a single simulation  
> parallel -j 2 slim -d s={3} -d mu={1} -d N=10000 '-d h="{2}"' -d Rep=1 -d r=1.0e-07 ScriptNeutralInversion_DefineInitialStateAutosome_AllDistrib.slim ::: 1.0e-08 ::: '"Distrib2"' ::: -0.0001 -0.001

*mu* define the mutation rate, *r* the recombination rate, *s* the selection coefficient (fixed or the mean of the gamma distribution), *h* the dominance coefficient, *N* the population size

Three options are possible for h. If h is a float, it will be considered constant for all mutations. If h is "Distrib1", mutation will have their selective coefficient randomly sampled among 0, 0.01, 0.1, 0.2, 0.3,0.4 and 0.5 with equal probabilities. If h is "Distrib2", mutations will have their selective coefficient randomly sampled among  0.1, 0.2, 0.3,0.4 and 0.5 with equal probabilities.

This writes the initial populations in a ../InitialState directory (must be created before!)

Other simulation 'burn-in' scripts are used for other simulations:
 - ScriptNeutralInversion_DefineInitialState_HaploDiplo_SingleChrom.slim #Simulation with HaploDiplontic life cycle ; From Jay et al. 2022
 - ScriptNeutralInversion_DefineInitialState_TwoChromFus.slim #Simulation of chromosome fusion (Figure S18); From Jay et al. 2022
 - ScriptNeutralInversion_DefineInitialState_XY_Variable_sh.slim #Simulation with mutation of variable h and s (drawn from distribution); From Jay et al. 2022
 - ScriptNeutralInversion_DefineInitialState.slim #Simulation of populations with non-XY system but with two mating types (all individuals are always heterozygous) (Figure S14)

## Introduce inversions
Then, we introduce in these initial populations 1000's of inversions or other recombination modifiers. To do so, use the snakemake script:
 - SnakeMakeSlimSimulationV5 --jobs 20 --configfile ParamFile_Distrib.yaml --config outputDir="Output/" initDir="InitialState/
 The "ParamFile_Distrib.yaml" file indicates the range of parameters to be simulated (see ParamFile_Distrib.yaml for an example) and the number of inversion to be simulated.
 This script run "Script_IntroInv_XY_AllChrom_AllH_SingleRepV2.slim", which allow to perform most simulations.

Other (older) script from Jay et al. 2022 allow to perform others type of simulations
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_SingleChrom_HaploDiplo.slim #Simulation with HaploDiplontic life cycle (Figures S22-23); From Jay et al. 2022
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_Variablesh.slim #Simulation with mutation of variable h and s (drawn from distribution, Figure S17); From Jay et al. 2022
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_ChromFus.slim #Simulation of chromosome fusion (Figure S18); From Jay et al. 2022
 - Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_HomoHeteroRecombMod.slim #Simulation of populations with non-XY system but with two mating types (all individuals are always heterozygous) and the introduction of not inversion but recombination modifier suppresssing recombination whatever their genotype (Figure S14)

Basic usage of those "IntroduceInversionFromInitStat" scripts from Jay et al. 2022:
 > slim -d N=1000 -d mu=1e-08 -d h=0.01 -d s=-0.01 -d r=1e-5 -d rep=1 -d start=1 -d end=2000000 -d "Init='../InitialState/slim_g15000_N=1000_r=1e-05_u=1e-08_s=-0.01_h=0.01_1635828730194.txt'" -d "SexChrom='Y'" Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY.slim #Introduction of a single 2000000 bp inversion on the Y chromosome in a initial population. Rep is only used to annotate the simulation (simunation number "Rep"). 

The position of the inversion is defined with *start* and *end*. To mimick chromosome fusion, just consider start and end spanning the gap between the two chromosome, *i.e.* the position 10Mb (*e.g. start=7000000, end=11000000)

Those scripts can be run with a modified snakemake file (not provided) or using, for instance, GNU parallel. 

This scripts produces output in ../Output (must be created) such as IntroInv_XYsyst_5M_SexChrom_N10000_r1.0e-07_u1.0e-09_s-0.0001_h0_sInv0_invSize1000000_All50000Rep_Stat.txt #file containing the trajectory of 50,000 1Mb sex-linked inversions with N=10000, s=-0.0001, h=0, etc # Result of the  Script_IntroInv_XY_AllChrom_AllH_SingleRepV2.slim script. Other script produce output with slightly different names.
For plotting, concatenate the output
 > cat cat IntroInv_XYsyst_5M_*N1* >> IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat_Full.txt 

 # Analyses of the formation of sex chromosomes (Figures 5-6 and S24):
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
