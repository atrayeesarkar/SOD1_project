#!/bin/bash

WT_c="/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_0/157_2"
WT_u='/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Uncharged/phic_0/157_2'
WT_m='/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_0/157_2'

mc_c='/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/157_2/rg_WT_c.xvg'
mc_u='/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_2/rg_WT_u.xvg'
mc_m='/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/157_2/rg_WT_m.xvg'

# cd $WT_c

# echo -e "0\n0" |gmx gyrate -f traj_comp_pbcmolcenter_WT_c.xtc -s topol.tpr -o $mc_c 

cd $WT_u

echo -e "0\n0" |gmx gyrate -f traj_comp_pbcmolcenter_WT_u.xtc -s topol.tpr -o $mc_u 

# cd $WT_m

# echo -e "0\n0" |gmx gyrate -f traj_comp_pbcmolcenter_WT_m.xtc -s topol.tpr -o $mc_m 