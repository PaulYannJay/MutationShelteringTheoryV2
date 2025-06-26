#!/bin/bash

slimFc(){ #Define a function
Init=`ls ../../InitialState/slim_g15000_MidSDR_10MChrom_XY_N=$1_r=1e-06_u=$2_s=$4_h=$3_*[0-9].txt` #Grep the Initial population file
position=$6 #Mid position of the inversion
size=$7 #Size of the inversion
end=$((position+size/2+1))
start=$((position-size/2+1))
~/Software/SLiM/build/slim -d N=$1 -d mu=$2 -d h=$3 -d s=$4 -d r=1e-6 -d rep=$5 -d start=$start -d end=$end -d "Init='$Init'" -d "SexChrom='$8'" Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_ModForSlim4.3.slim  #Launch script
}



