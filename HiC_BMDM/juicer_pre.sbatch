#!/bin/bash

# Script for generating a binned, normalized HiC file using juicer
# 'mm10_Arima.txt' generated from 'run_generate_site_positions.sbatch'
# 'mm10.chrom.sizes' supplied in repos

# Example usage:
# sbatch --array 0-2 juicer_pre_KR.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J pre
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Define query files
queries=($(ls *.sites.filtered.fragment.pairs | xargs -n 1 basename | sed 's/.sites.filtered.fragment.pairs//g' | uniq))

# Define temporary directory
tmpDir=juicer_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Generate .HiC file with high-quality pairs (MAPQ>30)
echo $(date +"[%b %d %H:%M:%S] Generating a binned, normalized HiC file using juicer (MAPQ>30)...")

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

java \
-Xmx64G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar juicer_tools.jar pre \
-f mm10_Arima.txt \
-r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 \
-k KR \
-q 30 \
-t ${tmpDir} \
-j ${numThreads} \
${queries[$SLURM_ARRAY_TASK_ID]}.sites.filtered.fragment.pairs \
${queries[$SLURM_ARRAY_TASK_ID]}.q30.KR.hic \
mm10.chrom.sizes

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -rf ${tmpDir}/

echo $(date +"[%b %d %H:%M:%S] Done")