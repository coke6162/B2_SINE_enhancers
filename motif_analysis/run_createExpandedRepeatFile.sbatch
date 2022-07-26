#!/bin/bash

# Script for "expanding" B2 coordinates relative to consensus according to RepeatMasker output
# repeats_turnRMFile_to_UCSC.py and createExpandedRepeatFile.py utility scripts available through repos
# Expanded bed files are available through repos
# RepeatMasker output file is availble here: https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm406-dfam2.0/mm10.fa.out.gz

# Example usage:
# sbatch run_createExpandedRepeatFile.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J expanded_repeats
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module purge
module load python/2.7.14

# Convert RepeatMasker out to UCSC tabular format
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Convert RepeatMasker output to UCSC tabular...")

repeats_turnRMFile_to_UCSC.py -i mm10_dfam.out > mm10_dfam.out.ucsc.txt

# Make "expanded" repeat bed files
echo $(date +"[%b %d %H:%M:%S] Expand repeats...")

createExpandedRepeatFile.py -i mm10_dfam.out.ucsc.txt | grep "B2_Mm2" | bedtools sort -i - > B2_Mm2_expanded.bed
createExpandedRepeatFile.py -i mm10_dfam.out.ucsc.txt | grep "B2_Mm1a" | bedtools sort -i - > B2_Mm1a_expanded.bed
createExpandedRepeatFile.py -i mm10_dfam.out.ucsc.txt | grep "B2_Mm1t" | bedtools sort -i - > B2_Mm1t_expanded.bed

echo $(date +"[%b %d %H:%M:%S] Done")