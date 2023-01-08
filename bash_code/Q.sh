#!/bin/bash

WT_c="/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_0/170_3"
WT_u='/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Uncharged/phic_0/170_3'
WT_m='/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_0/170_3'
# # traj -pbc no jump 
# nj_c='/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/157_2/traj_comp_pbcnojump_WT_c.xtc'
# nj_u='/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_2/traj_comp_pbcnojump_WT_u.xtc'
# nj_m='/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/157_2/traj_comp_pbcnojump_WT_m.xtc'
# traj -pbc mol center
mc_c='/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/170_3/traj_comp_pbcmolcenter_WT_c.xtc'
mc_u='/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/170_3/traj_comp_pbcmolcenter_WT_u.xtc'
mc_m='/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/170_3/traj_comp_pbcmolcenter_WT_m.xtc'

# # conf.gro -pbc no jump 
# cnj_c='/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/157_2/conf_pbcnojump_WT_c.gro'
# cnj_u='/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_2/conf_pbcnojump_WT_u.gro'
# cnj_m='/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/157_2/conf_pbcnojump_WT_m.gro'
# # conf.gro -pbc mol center
# cmc_c='/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/157_2/conf_pbcmolcenter_WT_c.gro'
# cmc_u='/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_2/conf_pbcmolcenter_WT_u.gro'
# cmc_m='/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/157_2/conf_pbcmolcenter_WT_m.gro'


cd $WT_c

#  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $nj_c -pbc nojump    
 echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_c -pbc mol -center -ur compact 
#  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $cnj_c -pbc nojump -e 0    
#  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $cmc_c -pbc mol -center -ur compact -e 0 
cd $WT_u

  # echo -e "0\n0" |gmx trjconv -f traj_comp.xtc -s topol.tpr -o $nj_u -pbc nojump  
  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_u -pbc mol -center -ur compact  
#   echo -e "0\n0" |gmx trjconv -f traj_comp.xtc -s topol.tpr -o $cnj_u -pbc nojump  -e 0
#   echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $cmc_u -pbc mol -center -ur compact  -e 0
cd $WT_m

#  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $nj_m -pbc nojump    
 echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_m -pbc mol -center -ur compact  
#  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $cnj_m -pbc nojump -e 0    
#  echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $cmc_m -pbc mol -center -ur compact  -e 0
