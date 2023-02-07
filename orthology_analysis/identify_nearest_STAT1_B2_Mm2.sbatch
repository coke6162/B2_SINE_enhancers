#!/bin/bash

# Script for running identifying the nearest STAT1-bound B2_Mm2 element relative to ISGs
# 'B2_Mm2_picc_IFNG_2h_STAT1.bed' and 'gencode.vM18.annotation.tss.bed' provided in repos

# Example usage:
# sbatch identify_nearest_STAT1_B2_Mm2.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J B2_proximity
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load module
module load bedtools/2.28.0

# Generate union ISG tss bed file
pwd; hostname; date

echo "Generating union ISG tss bed file..."

awk 'FNR==NR {a[$4]=$0; next} $1 in a {print a[$0]}' gencode.vM18.annotation.tss.bed \
picc_qiao_mouse_padj_log2FC0_up_IDs_union.txt \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss.bed

# Report distance to nearest STAT1-bound B2_Mm2 element
echo $(date +"[%b %d %H:%M:%S] For each union ISG, find nearest STAT1-bound B2_Mm2 element...")

bedtools sort -i picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss.bed \
| bedtools closest -d -t first -a - -b B2_Mm2_picc_IFNG_2h_STAT1.bed \
| sort -k13,13rg \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss_nearest_B2_Mm2_picc_IFNG_2h_STAT1.bed

# Filter for gene ID, distance to B2, and coordinates to nearest B2
echo $(date +"[%b %d %H:%M:%S] Condensing file...")

awk '{print $4 "\t" $13 "\t" "B2_Mm2_"$7":"$8"-"$9}' picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss_nearest_B2_Mm2_picc_IFNG_2h_STAT1.bed \
| sort -k1,1 \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss_nearest_B2_Mm2_picc_IFNG_2h_STAT1_condensed.txt

# Append to summarized binary matrix
echo $(date +"[%b %d %H:%M:%S] Appending to sumarized binary matrix...")

paste picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_summarized.txt \
picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss_nearest_B2_Mm2_picc_IFNG_2h_STAT1_condensed.txt \
| cut -f1-2,4-5 \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_nearest_B2_Mm2_picc_IFNG_2h_STAT1.txt

# Add gene name
echo $(date +"[%b %d %H:%M:%S] Adding gene name...")

awk 'NR==FNR{a[$1]=$0; next;} $1 in a {print a[$1] "\t" $2}' \
picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_nearest_B2_Mm2_picc_IFNG_2h_STAT1.txt \
picc_qiao_mouse_padj_log2FC0_up_IDs_names_union_binary.txt \
| awk '{OFS="\t"; print $1,$5,$2,$4,$3}' > tmp1.txt

# Add gene TSS coorinates
echo $(date +"[%b %d %H:%M:%S] Adding gene TSS coordinates...")

awk 'NR==FNR{a[$1]=$0; next;} $4 in a {print a[$4] "\t" $1":"$2"-"$3}' \
tmp1.txt picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss.bed \
| awk '{OFS="\t"; print $1,$2,$6,$3,$4,$5}' \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_nearest_B2_Mm2_picc_IFNG_2h_STAT1_summarized.txt

# Remove temporary files
echo $(date +"[%b %d %H:%M:%S] Removing temporary files...")

rm tmp*txt

echo $(date +"[%b %d %H:%M:%S] Done")