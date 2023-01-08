#!/bin/bash
id = (120,140,150,1551,156,157,158,159,160,170,180,190,200,210,220,230,240,250,260,280) 

for i in ${id[@]}; do
	scp asarkar4@uhpc.rcdc.uh.edu:/uhpc/cheung/taburt/SOD_project/METADATA/G41D_charged/phic0/$i_1/Etot_T.dat.gz /Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/i_1.txt.Z
done

