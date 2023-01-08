#!/bin/bash

#SBATCH -n 2
#SBATCH -N 1
#SBATCH -t 0-01:10
#SBATCH --mem-per-cpu=1GB
#SBATCH -o format%A.out
#SBATCH -e format%A.err

module purge
module load intel
./format.sh


