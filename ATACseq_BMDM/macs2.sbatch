#!/bin/bash

# Script for calling peaks with macs2

# Example usage:
# genomeSize=2.73e9 \
# sbatch --array=0-8 macs2.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=12:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J macs2
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Module load macs2
module load python/2.7.14/MACS/2.1.1

# Define query files
queries=($(ls *.sorted.dedup.shifted.bam | xargs -n 1 basename))

# Call peaks with macs2
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Calling peaks with macs2.."

macs2 callpeak \
--treatment ${queries[$SLURM_ARRAY_TASK_ID]} \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.shifted.bam} \
--outdir ./ \
--gsize 2.73e9 \
--SPMR -B --keep-dup all --format BAMPE --call-summits

echo $(date +"[%b %d %H:%M:%S] Done")