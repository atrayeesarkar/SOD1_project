#!/bin/bash

input="/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_20/160/160_1"

awk  'BEGIN{b=0; t=1}{a+=1; if($1=="ATOM")b =a; if($1=="END")t+=1; if((t>10000)&& (b==a)){fn="/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_20/160/J_textfiles/equil_"$6".txt"; print t"\t"$6"\t"$7"\t"$8"\t"$9 >>fn }}'$input
