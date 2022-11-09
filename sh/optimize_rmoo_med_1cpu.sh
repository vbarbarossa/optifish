#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=15G
#SBATCH -p "cpu-medium"
#SBATCH --output=sh/optimize_rmoo_med_1cpu_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/R-optifish2

Rscript R/optimize_rmoo.R
