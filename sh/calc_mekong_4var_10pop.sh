#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --array 1-10
#SBATCH -t 7-00:00:00
#SBATCH -p "cpu-long"
#SBATCH --output=sh/calc_mekong_4var.out
#SBATCH --mail-type="ALL"
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/r4

Rscript R/optimize_4var_10pop.R
