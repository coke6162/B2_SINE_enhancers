#!/bin/bash

## Script for running deeptools bamCoverage
## Example usage:
## inDir=path outDir=path sbatch --array 0-21 bamToBw.q

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=6:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bamToBw
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

module load singularity
binSize=1
numCPU=8

# Set query files
queries=($(ls $inDir/*.uniq.bam | xargs -n 1 basename))

# SRR8867628
singularity exec --bind /Shares/CL_Shared /scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img bamCoverage -b $inDir/${queries[$SLURM_ARRAY_TASK_ID]} -o $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.uniq.bam}_fwd.bw --binSize=$binSize -p $numCPU --normalizeUsing CPM --filterRNAstrand forward --ignoreForNormalization chrX chrM
singularity exec --bind /Shares/CL_Shared /scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img bamCoverage -b $inDir/${queries[$SLURM_ARRAY_TASK_ID]} -o $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.uniq.bam}_rev.bw --binSize=$binSize -p $numCPU --normalizeUsing CPM --filterRNAstrand reverse --ignoreForNormalization chrX chrM
