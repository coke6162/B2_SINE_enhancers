#!/bin/bash

# Script for merging sample WT & KO transcriptomes
# gtf_list.txt contains names of all sample GTFs (see repos for example)

# Example usage:
# sbatch stringtie_merge.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J stringtie_merge
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Define parameters
numThreads=8

# Load stringtie v1.3.3b
module load stringtie/1.3.3b

# Run stringtie
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Merging transcriptomes with Stringtie...")

stringtie --merge --rf -p ${numThreads} \
-o all_samples_merged.gtf \
-G gencode.vM18.annotation.gtf \
gtf_list.txt

echo $(date +"[%b %d %H:%M:%S] Done")