#!/bin/bash

## Script for running multiqc

## Example usage:
## inDir=. outDir=. sbatch multiqc.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J multiqc
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# load modules
module load singularity

# define key variables
multiqc=/scratch/Shares/public/singularity/multiqc-1.7.img

# Run MultiQC
pwd; hostname; date

echo "Singularity version: "$(singularity --version)
echo "MultiQC version: "$(singularity exec --bind /scratch ${multiqc} multiqc --version)
echo $(date +"[%b %d %H:%M:%S] Running MultiQC...")

singularity exec --bind /scratch ${multiqc} \
multiqc \
${inDir}/*.zip \
--outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")
