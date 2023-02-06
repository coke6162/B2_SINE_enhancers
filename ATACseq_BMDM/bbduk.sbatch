#!/bin/bash

# Script for running bbduk with paired-end data
# 'adapters.fa' supplied through github repos

# sbatch --array 0-8 bbduk.sbatch

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

# Define query samples
queries=($(ls *fastq.gz | xargs -n 1 basename | sed 's/_R1.fastq.gz//g' | sed 's/_R2.fastq.gz//g' | uniq))

# Run bbduk
pwd; hostname; date

echo "Processing sample: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running bbduk...")

bbduk.sh -Xmx8g in1=${queries[$SLURM_ARRAY_TASK_ID]}_R1.fastq.gz \
in2=${queries[$SLURM_ARRAY_TASK_ID]}_R2.fastq.gz \
out1=${queries[$SLURM_ARRAY_TASK_ID]}_trimmed_R1.fastq.gz \
out2=${queries[$SLURM_ARRAY_TASK_ID]}_trimmed_R2.fastq.gz \
ref=adapters.fa \
ktrim=r k=34 mink=11 hdist=1 tpe tbo \
qtrim=r trimq=10 \
t=${numThreads}

echo $(date +"[%b %d %H:%M:%S] Done")
