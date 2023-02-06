#!/bin/bash

# Script for adding juicer-style fragment information to pairs file using pairix fragment_4dnpairs.pl
# 'mm10_Arima.txt' generated from 'run_generate_site_positions.sbatch'

# Example usage:
# sbatch --array 0-2 run_addfrag2pairs.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J 4dnpairs
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Define query files
queries=($(ls *.filtered.pairs.gz | xargs -n 1 basename | sed 's/.filtered.pairs.gz//g' | uniq))

# Run fragment_4dnpairs.pl
pwd; hostname; date

echo "Running fragment_4dnpairs.pl..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

gunzip -c ${queries[$SLURM_ARRAY_TASK_ID]}.filtered.pairs.gz \
| fragment_4dnpairs.pl - \
${queries[$SLURM_ARRAY_TASK_ID]}.sites.filtered.fragment.pairs \
mm10_Arima.txt

echo $(date +"[%b %d %H:%M:%S] Done")