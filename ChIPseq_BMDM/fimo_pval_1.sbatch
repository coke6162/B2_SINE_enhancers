#!/bin/bash

# Script for identifying predicted TF motifs genome-wide
# PWM file in .meme format available through JASPAR
# For convenience PWM files are available through this repos
# "box_and_whisker/B2_all.fa" is output from "get_overlapping_B2_box_whisker.sbatch"

# Example usage:
# sbatch fimo_pval_1.sbatch

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

fimo --text --parse-genomic-coord --thresh 1 MA0137.3.meme box_and_whisker/B2_all.fa > box_and_whisker/GAS_FIMO_out_pval_1.txt

# Convert FIMO txt output to bed
pwd; hostname; date

awk '(NR>1){print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $8 "\t" $6}' box_and_whisker/GAS_FIMO_out_pval_1.txt | bedtools sort -i - > box_and_whisker/GAS_FIMO_out_pval_1.bed

echo $(date +"[%b %d %H:%M:%S] Done")