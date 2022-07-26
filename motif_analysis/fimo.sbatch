#!/bin/bash

# Script for identifying predicted TF motifs genome-wide
# PWM file in .meme format available through JASPAR
# For convenience PWM files are available through this repos

# Example usage:
# sbatch fimo.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J fimo
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Run FIMO
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Running FIMO...")

fimo --text --parse-genomic-coord --thresh 1e-4 MA0137.3.meme mm10.fa > GAS_fimo_out.txt
fimo --text --parse-genomic-coord --thresh 1e-4 MA0139.1.meme mm10.fa > CTCF_fimo_out.txt

echo $(date +"[%b %d %H:%M:%S] Done")