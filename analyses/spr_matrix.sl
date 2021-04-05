#!/bin/bash

#SBATCH --mail-user=kelly@ceremade.dauphine.fr
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o log/slurm-%j.out
#SBATCH -t 00:02:00

Rscript ../../src/spr_matrix.R $TARGET
