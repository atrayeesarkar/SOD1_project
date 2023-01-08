#!/bin/bash
## creating pbcmolcenter trajectories

#   for i in '120_1' '130_1' '140_1' '145_1' '150_1' '155_1' '158_1' '159_1' '160_1' '162_1' '163_1' '164_1' '165_1' '170_1' '180_1' '190_1' '200_1' '210_1' '220_1' 
# for i in '120_2' '140_2' '150_2' '155_2' '156_2' '157_2' '158_2' '159_2' '160_2' '170_2' '180_2' '190_2' '200_2' '210_2' '220_2' '230_2' '240_2' '250_2' '260_2' '280_2'
for i in '120_3' '140_3' '150_3' '155_3' '156_3' '157_3' '158_3' '159_3' '160_3' '180_3' 
 do
    WT_c="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/WT_charged/phic0/$i"
    WT_u="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/WT_uncharged/phic0/$i"
    WT_m="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/G41D_charged/phic0/$i"  
    mc_c="/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/$i/traj_comp_pbcmolcenter_WT_c.xtc"
    mc_u="/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/$i/traj_comp_pbcmolcenter_WT_u.xtc"
    mc_m="/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/$i/traj_comp_pbcmolcenter_WT_m.xtc"
    
   cd $WT_c
   # echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_c -pbc mol -center -ur compact  
   cp topol.tpr /Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/$i/
   cd $WT_u
   # echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_u -pbc mol -center -ur compact  
   cp topol.tpr /Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/$i/
   cd $WT_m
   # echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_m -pbc mol -center -ur compact  
   cp topol.tpr /Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/$i/
  
done

# WT_c="/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/WT_charged/phic20/161_1"
# WT_u='/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/WT_uncharged/phic20/161_1'
# WT_m='/Volumes/TAB_RESEARCH_CHEUNG/SOD_project/SIMULATIONS/G41D_charged/phic20/161_1'
# mc_c='/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_20/161_1/161_1.pdb'
# mc_u='/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_20/161_1/161_1.pdb'
# mc_m='/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_20/161_1/161_1.pdb'


# cd $WT_c
# echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_c 

# cd $WT_u
# echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_u 

# cd $WT_m
# echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -o $mc_m -pbc mol -center -ur compact  


# input="Data/WT_Uncharged/phic_20/161_1/traj_comp_pbcmolcenter_WT_u.xtc"
# $SKIPFRAMES=1;
# GMXPATH="/usr/local/Cellar/gromacs/2021.2/bin";
# print "What is the name of the tpr file?\n";
# $TPR=<STDIN>;
# chomp($TPR);
# # go through all of the xtc files
# print "What  xtc (or trr) file would you like to analyze?\n";
# $XTCfile=<STDIN>;
# chomp($XTCfile);
# # 
# `echo 0 | $GMXPATH/gmx trjconv  -skip $SKIPFRAMES -s $TPR -o tmp.pdb -f  $XTCfile`;
# input = "tmp.pdb";
# awk  'BEGIN{b=0;}{a+=1; if($2=="frame")t = $3+1; if($3==0)b=a; if((t>10000)&& (b==a)){fn="Data/WT_Uncharged/phic_20/157/J_textfiles/equil_"$3".txt"; print t"\t"$3"\t"$7"\t"$9"\t"$11 >>fn }}' $input #delete pre existing files of same names as the output files. ">>" append to the already existing files. 

# 'rm tmp.pdb`;

###### create directries
# cd ~/Desktop/SOD_project/Data/WT_Charged/phic_0/
# # mkdir { 120_1 130_1 140_1 150_1 155_1 158_1 159_1 160_1 162_1 163_1 164_1 165_1 170_1 180_1 190_1 200_1 210_1 220_1 }
# mkdir { 120_3 140_3 150_3 155_3 156_3 157_3 158_3 159_3 160_3 170_3 180_3 190_3 200_3 210_3 220_3 230_3 240_3 250_3 260_3 280_3 }

# cd ~/Desktop/SOD_project/Data/WT_Uncharged/phic_0/
# # mkdir { 120_1 130_1 140_1 150_1 155_1 158_1 159_1 160_1 162_1 163_1 164_1 165_1 170_1 180_1 190_1 200_1 210_1 220_1 }
# mkdir { 120_3 140_3 150_3 155_3 156_3 157_3 158_3 159_3 160_3 170_3 180_3 190_3 200_3 210_3 220_3 230_3 240_3 250_3 260_3 280_3 }

# cd ~/Desktop/SOD_project/Data/G41D_mutant/phic_0/
# # mkdir { 120_1 130_1 140_1 150_1 155_1 158_1 159_1 160_1 162_1 163_1 164_1 165_1 170_1 180_1 190_1 200_1 210_1 220_1 }
# mkdir { 120_3 140_3 150_3 155_3 156_3 157_3 158_3 159_3 160_3 170_3 180_3 190_3 200_3 210_3 220_3 230_3 240_3 250_3 260_3 280_3 }
