#!/bin/bash

# Script for running identifying B2_Mm2 elements predicted to interact with an ISG TSS
# 'gencode.vM18.annotation.tss.bed' & 'B2_Mm2.bed' provided in repos
# 'EnhancerPredictionsAllPutative.bed' generated through ABC_BMDM repos analysis

# Example usage:
# sbatch identify_interacting_ABC_B2_Mm2.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J B2_ABC
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load module
module load bedtools/2.28.0

# Extract all putative ABC enhancer/gene contacts with an ABC score >= 0.01
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Extracting enhancer/gene contacts (ABC >= 0.01)...")

awk '{if ($5 >= 0.01) print}' EnhancerPredictionsAllPutative.bed > EnhancerPredictionsAllPutative_0.01.bed

# Subset for ABC enhancer/gene contacts that contact an ISG TSS
echo $(date +"[%b %d %H:%M:%S] For each union ISG, find nearest STAT1-bound B2_Mm2 element...")

cut -f2 picc_qiao_mouse_padj_log2FC0_up_IDs_names_union.txt \
| awk 'FNR==NR{a[$4]=(a[$4]?a[$4] ORS:"")$0;next} ($1 in a){print a[$0]}' EnhancerPredictionsAllPutative_0.01.bed - \
> EnhancerPredictionsAllPutative_0.01_picc_qiao_mouse_padj_log2FC0_up.bed

# Intersecting with B2_Mm2 (requires the entire B2_Mm2 element to overlap)
echo $(date +"[%b %d %H:%M:%S] Intersecting with B2_Mm2...")

bedtools sort -i EnhancerPredictionsAllPutative_0.01_picc_qiao_mouse_padj_log2FC0_up.bed \
| bedtools intersect -a - -b B2_Mm2.bed -F 1 -wa -wb \
> EnhancerPredictionsAllPutative_0.01_picc_qiao_mouse_padj_log2FC0_up_B2_Mm2.bed

# Summarize into gene name, B2 coordinates, and ABC score
echo $(date +"[%b %d %H:%M:%S] Summarizing...")

awk '{if ($5 != "NaN") print $4 "\t" $10 "_" $7 ":" $8 "-" $9 "\t" $5}' \
EnhancerPredictionsAllPutative_0.01_picc_qiao_mouse_padj_log2FC0_up_B2_Mm2.bed \
> tmp1.txt

# Retain only top enhancer-gene contact
echo $(date +"[%b %d %H:%M:%S] Retaining only top enhancer-gene contact...")

sort -k3,3rn tmp1.txt | sort -k1,1 -u > tmp2.txt

# Convert to binary matrix (B2 genes only)
echo $(date +"[%b %d %H:%M:%S] Converting to binary column (B2 genes only)...")

awk 'NR==FNR{a[$2]=$0; next;} $1 in a {print a[$1] "\t" $2 "\t" $3}' \
picc_qiao_mouse_padj_log2FC0_up_IDs_names_union_binary.txt tmp2.txt \
> tmp3.txt

# Get binary matrix for genes without B2
echo $(date +"[%b %d %H:%M:%S] Converting to binary column (other genes)...")

cut -f1 tmp2.txt \
| grep -wvf - picc_qiao_mouse_padj_log2FC0_up_IDs_names_union_binary.txt \
| awk '{print $0 "\t" "NA" "\t" "NA"}' \
> tmp4.txt

# Concatenate
echo $(date +"[%b %d %H:%M:%S] Concatenating to form complete matrix...")

cat tmp3.txt tmp4.txt \
sort -k1,1 \
> tmp5.txt

# Add gene TSS coordinates
echo $(date +"[%b %d %H:%M:%S] Adding gene TSS coordinates...")

awk 'NR==FNR{a[$1]=$0; next;} $4 in a {print a[$4] "\t" $1":"$2"-"$3}' tmp5.txt picc_qiao_mouse_padj_log2FC0_up_IDs_union_tss.bed \
| awk '{OFS="\t"; print $1,$2,$6,$3,$4,$5}' \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_B2_Mm2_ABC_0.01_summarized.txt

# Remove temporary files
echo $(date +"[%b %d %H:%M:%S] Removing temporary files...")

rm tmp*txt

echo $(date +"[%b %d %H:%M:%S] Done")