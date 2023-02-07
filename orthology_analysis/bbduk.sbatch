#!/bin/bash

# Script for running bbduk with single-end unstranded data
# 'adapters.fa' supplied through github repos

# Example usage:
# sbatch --array 0-3 bbduk.sbatch.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J bbduk
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=1

# Load module
module load bbmap/38.05

# Define query files
queries=($(ls *.fastq.gz | sed 's/.fastq.gz//g' | uniq))

# Define adapter file
adapterFile=/opt/bbmap/38.05/resources/adapters.fa

# Run bbduk
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running bbduk...")

bbduk.sh -Xmx4g \
in=${queries[$SLURM_ARRAY_TASK_ID]}.fastq.gz \
out=${queries[$SLURM_ARRAY_TASK_ID]}_trimmed.fastq.gz \
ref=adapters.fa \
ktrim=r k=34 mink=11 hdist=1 tpe tbo \
qtrim=r trimq=10 \
t=${numThreads}

echo $(date +"[%b %d %H:%M:%S] Done")