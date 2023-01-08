#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2-20:30
#SBATCH -o Jaccard%A.out
#SBATCH -e Jaccard%A.err
#SBATCH --mem 125GB
#SBATCH --cpus-per-task=12

module load intel
module load matlab

matlab -nodisplay -nosplash -r P_inter_unfolding_folding_crowded_v3


