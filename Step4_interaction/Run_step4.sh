#!/bin/bash
#SBATCH --job-name=EF_psy_corr_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=70
#SBATCH --mem-per-cpu=20G 
#SBATCH -p q_fat_c
#SBATCH -q high_c

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

module load R/4.2.2
Rscript ./Step3_corr_slope.R