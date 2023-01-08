#!/bin/bash
#  for i in '120_1' '130_1' '140_1' '145_1' '150_1' '155_1' '158_1' '159_1' '160_1' '161_1' '162_1' '163_1' '164_1' '165_1' '170_1' '180_1' '190_1' '200_1' '210_1' '220_1' 
# for i in '120_2' '140_2' '150_2' '155_2' '156_2' '157_2' '158_2' '159_2' '160_2' '170_2' '180_2' '190_2' '200_2' '210_2' '220_2' '230_2' '240_2' '250_2' '260_2' '280_2'
 for i in '120_3' '140_3' '150_3' '155_3' '156_3' '157_3' '158_3' '159_3' '160_3' '180_3'
 do
    WT_c="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/WT_charged/phic0/$i"
    WT_u="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/WT_uncharged/phic0/$i"
    WT_m="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/G41D_charged/phic0/$i"  
    mc_c="/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/$i/energy.xvg"
    mc_u="/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/$i/energy.xvg"
    mc_m="/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/$i/energy.xvg"
    
   cd $WT_c
   echo -e "10\n0" | gmx energy -f ener.edr -s topol.tpr -o $mc_c
   input=$mc_c
   output="/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/$i/energy.txt"
   awk  'BEGIN{}{if(($1!="@")&&($1!="@TYPE")&&($1!="#")){print $1"\t"$2}}' $input>$output
	
   cd $WT_u
   echo -e "10\n0" | gmx energy -f ener.edr -s topol.tpr -o $mc_u
   input=$mc_u
   output="/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/$i/energy.txt"
   awk  'BEGIN{}{if(($1!="@")&&($1!="@TYPE")&&($1!="#")){print $1"\t"$2}}' $input>$output
	
   cd $WT_m
   echo -e "10\n0" | gmx energy -f ener.edr -s topol.tpr -o $mc_m 
   input=$mc_m
   output="/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/$i/energy.txt"
   awk  'BEGIN{}{if(($1!="@")&&($1!="@TYPE")&&($1!="#")){print $1"\t"$2}}' $input>$output
	

done