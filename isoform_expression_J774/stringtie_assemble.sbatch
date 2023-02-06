#!/bin/bash

# Script for assembling transcripts with Stringtie

# Example usage:
# sbatch --array 0-29 stringtie_assemble.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --time=1:00:00
#SBATCH --mem=2GB

# Job name and output
#SBATCH -J stringtie_assemble
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=2

# Load stringtie v1.3.3b
module load stringtie/1.3.3b

# Define query files
queries=($(ls *.uniq.bam | xargs -n 1 basename | sed 's/.uniq.bam//g'))

# Run stringtie
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}.uniq.bam
echo $(date +"[%b %d %H:%M:%S] Assembling transcripts with Stringtie...")

stringtie \
{queries[$SLURM_ARRAY_TASK_ID]}.uniq.bam \
-o ${queries[$SLURM_ARRAY_TASK_ID]}.gtf \
-G gencode.vM18.annotation.gtf \
--rf -j 5 -p ${numThreads} 

echo $(date +"[%b %d %H:%M:%S] Done")