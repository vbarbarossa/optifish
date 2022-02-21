#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH -p "cpu-long"
#SBATCH --output=sh/calc_mekong_4var_lev12.out
#SBATCH --mail-type="ALL"
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/r4

Rscript R/optimize_4var_lev12_merged.R
