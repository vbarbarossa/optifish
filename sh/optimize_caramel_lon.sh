#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 7-00:00:00
#SBATCH -p "cpu-long"
#SBATCH --output=sh/optimize_caramel_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=v.barbarossa@cml.leidenuniv.nl

module load Miniconda3
source activate ~/envs/R-optifish2

Rscript R/optimize_caramel_future_mitig.R
