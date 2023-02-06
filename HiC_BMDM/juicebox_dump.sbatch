#!/bin/bash

# Script for converting .hic to ABC-compatible format

# Example usage:
# sbatch --array 0-2 juicebox_dump.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J dump
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Load modules
module purge
module load gcc/7.1.0 bedtools python/3.7.4

# Define query files
queries=($(ls *.q30.KR.hic | xargs -n 1 basename | sed 's/.q30.KR.hic//g' | uniq))

# Define temporary directory
tmpDir=juicer_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Run juicebox dump
echo $(date +"[%b %d %H:%M:%S] Running juicebox dump...")

for i in {1..19} X Y
do
	mkdir -p ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/chr${i}
	mkdir -p ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/powerlaw
	java \
	-Xmx8G -Djava.io.tmpdir=${tmpDir} -XX:ParallelGCThreads=${numThreads} \
	-jar juicer_tools.jar dump \
	observed KR ${queries[$SLURM_ARRAY_TASK_ID]}.q30.KR.hic chr${i} chr${i} BP 5000 ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/chr${i}/chr${i}.KRobserved
	gzip ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/chr${i}/chr${i}.KRobserved
done

for i in {1..19} X Y
do
	mkdir -p ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/chr${i}
	mkdir -p ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/powerlaw
	java \
	-Xmx8G -Djava.io.tmpdir=${tmpDir} -XX:ParallelGCThreads=${numThreads} \
	-jar juicer_tools.jar dump \
	norm KR ${queries[$SLURM_ARRAY_TASK_ID]}.q30.KR.hic chr${i} BP 5000 ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/chr${i}/chr${i}.KRnorm
	gzip ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/chr${i}/chr${i}.KRnorm
done

# Fit HiC data to powerlaw model and extract parameters
echo $(date +"[%b %d %H:%M:%S] Fitting HiC data to powerlaw model...")

python /Users/coke6162/programs/ABC-Enhancer-Gene-Prediction-0.2.2/src/compute_powerlaw_fit_from_hic.py \
--hicDir ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr \
--outDir ${queries[$SLURM_ARRAY_TASK_ID]}/by_chr/powerlaw \
--maxWindow 1000000 \
--minWindow 5000 \
--resolution 5000 \
--chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -rf ${tmpDir}/

echo $(date +"[%b %d %H:%M:%S] Done")