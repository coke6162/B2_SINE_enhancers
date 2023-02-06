#!/bin/bash

# Script for computing ABC scores to predict enhancers with Hi-C data

# Example usage:
# sbatch predict_enhancers.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=32G

# Job name and output
#SBATCH -J predict_enhancers
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module purge
module load gcc/7.1.0 bedtools/2.28.0 python/3.7.4

# Predict enhancers
echo $(date +"[%b %d %H:%M:%S] Predicting enhancers...")

python ABC-Enhancer-Gene-Prediction-0.2.2/src/predict.py \
--enhancers EnhancerList.txt \
--genes GeneList.txt \
--HiCdir plat_BMDM_IFNG_2h_HiC/by_chr/ \
--hic_type juicebox \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold 0.02 \
--outdir IFNG \
--make_all_putative

echo $(date +"[%b %d %H:%M:%S] Done")