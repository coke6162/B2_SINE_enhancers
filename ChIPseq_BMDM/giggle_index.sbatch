#!/bin/bash

# Script to prepare GIGGLE index using mm10 Dfam TE annotation
# gzip'd "mm10_dfam_filtered.bed" file available in "RNAseq_BMDM" repos

# sbatch giggle_index.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --propagate=NONE

# Job name and output
#SBATCH -J giggle_index
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module load samtools/1.10

# Make directories for GIGGLE index
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Making directories...")

mkdir -p giggle_index/repeats
mkdir -p giggle_index/sorted
mkdir -p giggle_index/indexed

# Remove simple and tRNA repeats from .bed
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Removing simple and tRNA repeats...")

grep -v ")n" mm10_dfam.bed \
| grep -E 'SINE|LINE|LTR|DNA' \
| grep -v "\?" \
> mm10_dfam_filtered.bed

# Create a file that contains all unique TE family names
echo $(date +"[%b %d %H:%M:%S] Create a txt file containing all unique TE family names from filtered bed...")

awk '{print $5}' mm10_dfam_filtered.bed \
| sort \
| uniq \
> mm10_dfam_filtered_names.txt

# Extract repeats, writing one bed file per family to the "giggle_index/repeats/" subdir
echo $(date +"[%b %d %H:%M:%S] Extracting repeats...")

awk 'NR==FNR{ids[$1]; next} $4 in ids{print > ("giggle_index/repeats/"$4".bed")}' mm10_dfam_filtered_names.txt mm10_dfam_filtered.bed

# Compress bed files
echo $(date +"[%b %d %H:%M:%S] Compressing repeat files...")

for i in giggle_index/repeats/*.bed
do
	bgzip ${i}
done

# Sort compressed bed files
echo $(date +"[%b %d %H:%M:%S] Sorting compressed repeat files...")

sort_bed.sh "giggle_index/repeats/*.bed.gz" giggle_index/sorted 4

# Create a giggle index
echo $(date +"[%b %d %H:%M:%S] Index sorted, compressed bed files...")

time giggle index -i "giggle_index/sorted/*.bed.gz" -o giggle_index/indexed -f -s

echo $(date +"[%b %d %H:%M:%S] Done")