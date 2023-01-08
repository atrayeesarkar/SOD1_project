#!/bin/bash

#input="/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157/157_1"
input="/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_3/157_3.pdb"

#for i in {1..20}
#do
#output="/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_20/160/J_textfiles/equil_$i.txt"
#awk -v r=$i 'BEGIN{b=0; t=1}{a+=1; if($1=="ATOM")b =a; if($1=="END")t+=1; if((t>=10000)&& (b==a) && ($6==r)) print t"\t"$6"\t"$7"\t"$8"\t"$9}' $input > $output &
#done

#wait

awk  'BEGIN{b=0; t=1}{a+=1; if($1=="ATOM")b =a; if($1=="TER")t+=1; if((t>10000)&& (b==a)){fn="/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_3/J_textfiles/equil_"$5".txt"; print t"\t"$5"\t"$6"\t"$7"\t"$8 >>fn }}' $input #delete pre existing files of same names as the output files. ">>" append to the already existing files. 

