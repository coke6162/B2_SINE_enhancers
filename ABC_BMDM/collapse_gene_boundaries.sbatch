#!/bin/bash

# Script for collapsing transcripts bed to gene boundary bed
# 'gencode.vM18.annotation.bed' supplied in repos

# Example usage:
# sbatch collapse_gene_boundaries.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J collapse
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module load bedtools/2.28.0

# Navigate to working directory
pwd; hostname; date

echo "Navigating to working directory..."

cd ${inDir}/

# Make temporary directory
echo $(date +"[%b %d %H:%M:%S] Making temporary directory...")

mkdir tmp/

# Collapse transcripts bed to gene boundary bed
echo $(date +"[%b %d %H:%M:%S] Collapsing transcripts bed to gene boundary bed...")

for i in $(awk '{print $4}' gencode.vM18.annotation.bed | sort | uniq)
do
	awk -v name="$i" '{if ($4 == name) print $0}' gencode.vM18.annotation.bed > tmp/${i}.bed.tmp
	for j in $(awk '{print $1}' tmp/${i}.bed.tmp  | uniq)
	do
		awk -v chr="$j" '{if ($1 == chr) print $0}' tmp/${i}.bed.tmp > tmp/${i}_${j}.bed.tmp
		for k in $(awk '{print $6}' tmp/${i}_${j}.bed.tmp)
		do
			awk -v strand="$k" '{if ($6 == strand) print $0}' tmp/${i}_${j}.bed.tmp > tmp/${i}_${j}_${k}.bed.tmp
			start=$(sort -k2,2n tmp/${i}_${j}_${k}.bed.tmp | awk '(NR==1){print $2}')
			end=$(sort -k3,3rn tmp/${i}_${j}_${k}.bed.tmp | awk '(NR==1){print $3}')
			awk -v start="$start" -v end="$end" '{print $1 "\t" start "\t" end "\t" $4 "\t" $5 "\t" $6}' tmp/${i}_${j}_${k}.bed.tmp
			rm tmp/${i}_${j}_${k}.bed.tmp
		done | sort -k2,2n -u
		rm tmp/${i}_${j}.bed.tmp
	done
	rm tmp/${i}.bed.tmp
done | bedtools sort -i - > gencode.vM18.annotation.boundaries.bed

# Retain only one entry per gene (otherwise ABC will error)
echo $(date +"[%b %d %H:%M:%S] Retaining only one entry per gene...")

sort -k4,4 -u gencode.vM18.annotation.boundaries.bed | bedtools sort -i - > gencode.vM18.annotation.boundaries.uniq.bed

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -r tmp/

echo $(date +"[%b %d %H:%M:%S] Done")