#!/bin/bash

# Script to generate restriction site file for Arima enzyme
# 'generate_site_positions.py' is available through Juicer
# For convenience have included 'generate_site_positions.py' in repos (from Juicer v1.6)
# 'mm10.fa' not supplied, see README for link to download

# Example usage:
# sbatch run_generate_site_positions.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J RE_site
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module purge
module load python/3.7.14

# Generate restriction site file
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Generating restriction site file...")

generate_site_positions.py Arima mm10.fa

echo $(date +"[%b %d %H:%M:%S] Done")