#!/bin/bash

# Script for generating custom TEtranscripts GTF
# STAT1-bound B2_Mm2 elements are annotated separately from unbound B2_Mm2
# B2_Mm2 status is defined using paired STAT1 ChIP-seq data

# This script takes RepeatMasker output generated using Dfam v2.0 annotation
# RepeatMasker output file is availble here: https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm406-dfam2.0/mm10.fa.out.gz
# gzip'd intermediate bed file ("mm10_dfam.bed") available on repos

# For convenience, ChIP-seq files are provided on repos
# See "ChIPseq_BMDM" repos for how to generate these files

# Example usage:
# sbatch TEtranscripts.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=2-0:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J generate_gtf
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Convert RepeatMasker output file to bed format
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Converting RepeatMasker out to bed...")

awk '(NR>3){OFS="\t"; if ($11 ~ "/") print $5,$6-1,$7,$9,$10,$11,$1; else print $5,$6,$7,$9,$10,$11"/"$11,$1}' mm10_dfam.out \
| tr '/' '\t' \
| awk '{OFS="\t"; if ($4 == "C") print $1,$2,$3,"-",$5,$7,$6,$8; else print $1,$2,$3,$4,$5,$7,$6,$8}' \
| bedtools sort -i - \
> mm10_dfam.bed

# Reannotate STAT1-bound B2_Mm2 according to Piccolo et al. 2017 ChIP-seq
echo $(date +"[%b %d %H:%M:%S] Reannotate B2_Mm2 (Piccolo)...")

bedtools intersect -a mm10_dfam.bed -b B2_Mm2_piccolo_IFNG_2h_STAT1_intersected.bed -wao -r -f 1 \
| awk '{OFS="\t"; if ($9 != "." && $5 == "B2_Mm2") print $1,$2,$3,$4,"STAT1_"$5,$6,$7,$8; else if ($9 == "." && $5 == "B2_Mm2") print $1,$2,$3,$4,"Unbound_"$5,$6,$7,$8; else if ($9 == "." && $5 != "B2_Mm2") print $1,$2,$3,$4,$5,$6,$7,$8}' \
> mm10_dfam_B2_Mm2_picc_IFNG_2h_STAT1_modified.bed

# Reannotate STAT1-bound B2_Mm2 according to Platanitis et al. 2019 ChIP-seq
echo $(date +"[%b %d %H:%M:%S] Reannotate B2_Mm2 (Platanitis)...")

bedtools intersect -a mm10_dfam.bed -b B2_Mm2_plat_IFNG_1.5h_STAT1_intersected.bed -wao -r -f 1 \
| awk '{OFS="\t"; if ($9 != "." && $5 == "B2_Mm2") print $1,$2,$3,$4,"STAT1_"$5,$6,$7,$8; else if ($9 == "." && $5 == "B2_Mm2") print $1,$2,$3,$4,"Unbound_"$5,$6,$7,$8; else if ($9 == "." && $5 != "B2_Mm2") print $1,$2,$3,$4,$5,$6,$7,$8}' \
> mm10_dfam_B2_Mm2_plat_IFNG_1.5h_STAT1_modified.bed

# Download utility script to convert bed to gtf
echo $(date +"[%b %d %H:%M:%S] Download makeTEgtf.pl utility script...")

wget http://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/makeTEgtf.pl.gz
gzip -d makeTEgtf.pl.gz
chmod 770 makeTEgtf.pl

# Generate custom GTF referencing Piccolo et al. 2017 STAT1-bound B2_Mm2
echo $(date +"[%b %d %H:%M:%S] Generate custom GTF (Piccolo)...")

./makeTEgtf.pl \
-c 1 -s 2 -e 3 -o 4 \
-n mm10_dfam -t 5 -f 6 -C 7 -S 8 \
mm10_dfam_B2_Mm2_picc_IFNG_2h_STAT1_modified.bed \
> mm10_dfam_B2_Mm2_picc_IFNG_2h_STAT1_modified.gtf

# Generate custom GTF referencing Platanitis et al. 2019 STAT1-bound B2_Mm2
echo $(date +"[%b %d %H:%M:%S] Generate custom GTF (Platanitis)...")

./makeTEgtf.pl \
-c 1 -s 2 -e 3 -o 4 \
-n mm10_dfam -t 5 -f 6 -C 7 -S 8 \
mm10_dfam_B2_Mm2_plat_IFNG_1.5h_STAT1_modified.bed \
> mm10_dfam_B2_Mm2_plat_IFNG_1.5h_STAT1_modified.gtf

echo $(date +"[%b %d %H:%M:%S] Done")