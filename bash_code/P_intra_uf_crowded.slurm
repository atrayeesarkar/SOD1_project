#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-20:30
#SBATCH -o Jaccard%A.out
#SBATCH -e Jaccard%A.err
#SBATCH --mem-per-cpu=5GB
#SBATCH --cpus-per-task=20 

module load intel
module load matlab

matlab -nodisplay -nosplash -r P_intra_uf_crowded
~                                  
