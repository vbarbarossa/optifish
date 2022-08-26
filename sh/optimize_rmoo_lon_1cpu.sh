#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -p "cpu-long"
#SBATCH --output=sh/optimize_rmoo_lon_1cpu_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/R-optifish2

Rscript R/optimize_rmoo.R
