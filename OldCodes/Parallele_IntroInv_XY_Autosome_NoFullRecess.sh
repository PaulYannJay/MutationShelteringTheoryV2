#!/bin/bash

slimFc(){ #Define a function
Init=`ls ../InitialState/slim_g200000_XYsyst_5M_Autosome_N$1_r1e-8_u1e-08_s$2_NoFullRecess_Rep1_*[0-9].txt`
echo $Init
position=$3 #Mid position of the inversion
size=$4 #Size of the inversion
sInv=$5
end=$((position+size/2+1))
start=$((position-size/2+1))
~/Software/SLiM/build/slim -d N=$1 -d s=$2 -d sInv=$5 -d Rep=$6 -d start=$start -d end=$end -d "Init='$Init'" Script_IntroInv_XY_Autosome_NoFullRecess.slim #Launch script

}



