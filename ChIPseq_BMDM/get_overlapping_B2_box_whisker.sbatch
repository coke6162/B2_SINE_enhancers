#!/bin/bash

# Script for identifying B2 bound by STAT1, CTCF only, both STAT1 and CTCF, or neither for box-and-whisker plot
# "mm10_dfam.bed" is available in gzip'd format in the "RNAseq_BMDM" repos

# Example usage:
# sbatch get_overlapping_B2_box_whisker.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J overlapping_B2
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load bedtools
module load bedtools/2.28.0

# Extract B2 coordinates
grep "B2_Mm2" mm10_dfam.bed > B2_Mm2.bed
grep "B2_Mm1a" mm10_dfam.bed > B2_Mm1a.bed
grep "B2_Mm1t" mm10_dfam.bed > B2_Mm1t.bed

# Set shuffle seed
# Seeding adopted from https://stackoverflow.com/a/41962458/7820599
# Option parsing adopted from https://stackoverflow.com/a/14203146
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

seed=1;

REST=""
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	-s)
	    seed="$2"
	    shift
	    shift
	    ;;
	*)    # unknown option
	    REST="$REST $1"
	    shift # past argument
	    ;;
    esac
done

# Make subdirectory
mkdir -p box_and_whisker/

# Identify top unique GAS motif (by p-value) overlapping IFNG STAT1-bound B2 elements
for i in B2_Mm2
do
	bedtools intersect -a ${i}.bed -b picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -u | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_"NR "\t" $5 "\t" $6}' > box_and_whisker/${i}_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed
done

# Identify UT CTCF-bound B2 elements
for i in B2_Mm2 B2_Mm1a B2_Mm1t
do
	bedtools intersect -a ${i}.bed -b gual_BMDM_shNT_UT_CTCF_intersected.bed -u > box_and_whisker/${i}_gual_BMDM_shNT_UT_CTCF_intersected.bed
done

# Identify top unique GAS motif (by p-value) overlapping random subset of 5000 control B2 elements that don't overlap CTCF or STAT1
for i in B2_Mm2 B2_Mm1a B2_Mm1t
do
	bedtools intersect -a ${i}.bed -b box_and_whisker/${i}*_intersected.bed -u | awk '{if ($1 != "chrY") print}' | shuf -n 5000 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_"NR "\t" $5 "\t" $6}' > box_and_whisker/${i}_control_no_STAT1_or_CTCF_shuf5000.bed
done

# Identify B2_Mm2 bound by CTCF but not STAT1
bedtools intersect -a box_and_whisker/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected.bed -b box_and_whisker/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -u | awk '{if ($1 != "chrY") print}' | shuf -n 5000 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_"NR "\t" $5 "\t" $6}' > box_and_whisker/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_no_STAT1.bed

# Concatenate all bed files ending in "_GAS.bed" and convert to fasta
cat B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed B2_Mm2_control_no_STAT1_or_CTCF_shuf5000.bed B2_Mm1a_control_no_STAT1_or_CTCF_shuf5000.bed B2_Mm1t_control_no_STAT1_or_CTCF_shuf5000.bed B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_no_STAT1.bed | bedtools sort -i - | bedtools getfasta -fi mm10.fa -bed - -fo box_and_whisker/B2_all.fa