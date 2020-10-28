#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=compile.txt
#SBATCH --mem=5Gb
#SBATCH --mail-user=m.albecker@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=00:30:00
#SBATCH --output=compile.output
#SBATCH --error=compile.error

Rscript --vanilla Result_Compile.R
