#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 6-00:00:00
#SBATCH -p "cpu-long"
#SBATCH --output=sh/calc_mekong.out
#SBATCH --mail-type="ALL"
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/r4

Rscript R/optimize.R
