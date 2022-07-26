#!/bin/bash

## Script for calling peaks with macs2 (no input ctrls)

## Example usage:
## inDir=path outDir=path genome=mm sbatch --array 0-21 MACS2.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=12:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J MACS2
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# define query files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))

# load modules
module load python/2.7.14/MACS/2.1.1

# run the thing
pwd; hostname; date

echo "macs2 version: "$(macs2 --version)
echo "Target file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting macs2...")

macs2 callpeak \
--format BAMPE \
--treatment ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam} \
-g ${genome} \
--SPMR -B \
--outdir ${outDir} \
--call-summits

echo $(date +"[%b %d %H:%M:%S] Done!")
