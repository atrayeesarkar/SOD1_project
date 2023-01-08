#!/bin/bash

input="single_CA.pdb"
output="single_CA.txt"
awk  'BEGIN{b=0;}{a+=1; if($1=="ATOM")b =a;  if((b==a)&&($3=="CA")){ print $5"\t"$6"\t"$7"\t"$8 }} ' $input > $output



# input="crowd_CA.pdb"
# output="crowd_CA_1.txt"
# awk  'BEGIN{b=0;}{a+=1; if($1=="ATOM")b =a; if((b==a)&&($2<=110)){ print $5"\t"$6"\t"$7"\t"$8 }} ' $input > $output
