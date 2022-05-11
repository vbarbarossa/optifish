#!/bin/bash
#SBATCH -N 1
#SBATCH --array 1-16
#SBATCH -n 24
#SBATCH -t 1-00:00:00
#SBATCH -p "cpu-medium"
#SBATCH --output=sh/calc_mekong_sens_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/r4

Rscript R/calc_sensitivity.R
