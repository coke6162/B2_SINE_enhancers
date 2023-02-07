#!/bin/bash

# Script for generating a binary matrix reporting ISG induction by species
# 'gencode.vM18.annotation.tss.bed' provided in repos

# Example usage:
# sbatch generate_binary_matrix.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J binary
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load module
module load bedtools/2.28.0

# Extract ISG IDs
pwd; hostname; date

echo "Extracting gene IDs for all human ISGs..."

awk 'BEGIN { OFS=FS="\t" } (NR>1) { sub("\\..*", "", $1); print $1}' qiao_mono_IFNG_24h_padj_log2FC0_up.txt \
> qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable.txt

# Extract one2one, confidence score=1 orthologs
echo $(date +"[%b %d %H:%M:%S] Filtering for one2one, confidence score=1 orthologs...")

awk '{if ($5 == "ortholog_one2one" && $6 == 1) print}' human_ensembl_v105_biomart_to_mouse_stable.txt \
> human_ensembl_v105_biomart_to_mouse_stable_one2one_confidence1.txt

# Convert to mouse IDs (one2one, confidence score=1 only)
echo $(date +"[%b %d %H:%M:%S] Convert human ISGs to mouse orthologous IDs...")

awk -F "\t" 'NR==FNR{a[$1]=$0; next} ($1 in a){print a[$1]}' \
human_ensembl_v105_biomart_to_mouse_stable_one2one_confidence1.txt \
qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable.txt \
| cut -f4,5 \
> qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable_to_mouse_one2one_confidence1.txt

# Convert txt to bed
echo $(date +"[%b %d %H:%M:%S] Converting ISG orthologs txt to TSS bed...")

awk '{print $3}' qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable_to_mouse_one2one_confidence1.txt \
| awk 'FNR==NR {a[$4]=$0; next} $1 in a {print a[$0]}' gencode.vM18.annotation.tss.bed - \
> qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable_to_mouse_one2one_confidence1_tss.bed

# Ensure mouse ISG IDs do not have version number
echo $(date +"[%b %d %H:%M:%S] Removing gene ID version from mouse ISG ID...")

awk 'BEGIN { OFS=FS="\t" } { sub("\\..*", "", $4); print }' picc_BMDM_IFNG_4h_padj_log2FC0_up_tss.bed \
> picc_BMDM_IFNG_4h_padj_log2FC0_up_tss_stable.bed

# Define union mouse & human-to-mouse ISG gene ID set
echo $(date +"[%b %d %H:%M:%S] Defining union mouse & human-to-mouse ISG gene ID set...")

for i in qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable_to_mouse_one2one_confidence1_tss.bed picc_BMDM_IFNG_4h_padj_log2FC0_up_tss_stable.bed
do 
	cut -f4 ${i} 
done | sort | uniq > picc_qiao_mouse_padj_log2FC0_up_IDs_union.txt

# Add gene names to union list
echo $(date +"[%b %d %H:%M:%S] Adding gene names...")

awk -F "\t" 'NR==FNR{a[$4]=$4 "\t" $5; next} ($1 in a){print a[$0]}' \
gencode.vM18.annotation.tss.bed \
picc_qiao_mouse_padj_log2FC0_up_IDs_union.txt \
> picc_qiao_mouse_padj_log2FC0_up_IDs_names_union.txt

# Define which human-to-mouse ISGs are present in union
echo $(date +"[%b %d %H:%M:%S] Defining which human-to-mouse ISGs are present in union set...")

awk '{print $3}' qiao_mono_IFNG_24h_padj_log2FC0_up_IDs_stable_to_mouse_one2one_confidence1.txt \
| awk 'NR==FNR{a[$1]=$1; next;} $1 in a {delete a[$1]; print $1 "\t" 1;} END {for (key in a){print a[key] "\t" 0} }' \
picc_qiao_mouse_padj_log2FC0_up_IDs_union.txt - \
| sort -k1,1rn \
> human.txt.tmp

# Define which mouse ISGs are present in union
echo $(date +"[%b %d %H:%M:%S] Defining which mouse ISGs are present in union set...")

awk '{print $4}' picc_BMDM_IFNG_4h_padj_log2FC0_up_tss_stable.bed \
| awk 'NR==FNR{a[$1]=$1; next;} $1 in a {delete a[$1]; print $1 "\t" 1;} END {for (key in a){print a[key] "\t" 0} }' \
picc_qiao_mouse_padj_log2FC0_up_IDs_union.txt - \
| sort -k1,1rn \
> mouse.txt.tmp

# Generate complete binary matrix (col1=ID; col2=mouse binary; col3=human binary)
echo $(date +"[%b %d %H:%M:%S] Generating complete binary matrix...")

paste mouse.txt.tmp human.txt.tmp \
| cut -f1,2,4 \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary.txt

# Summarize to two columns (col1=ID; col2=species specificity)
# "mouse"=not induced according to human-to-mouse liftover
# "human"=not induced in mouse
# "shared"=induced in mouse and human-to-mouse
echo $(date +"[%b %d %H:%M:%S] Summarizing matrix...")

awk '{if ($2 == 1 && $3 == 0) print $1 "\t" "mouse"; else if($2 == 0 && $3 == 1) print $1 "\t" "human"; else print $1 "\t" "shared"}' picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary.txt \
> picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_summarized.txt

# Append gene names
echo $(date +"[%b %d %H:%M:%S] Appending gene names...")

paste picc_qiao_mouse_padj_log2FC0_up_IDs_names_union.txt picc_qiao_mouse_padj_log2FC0_up_IDs_union_binary_summarized.txt \
| cut -f1,2,4 \
> picc_qiao_mouse_padj_log2FC0_up_IDs_names_union_binary_summarized.txt

# Remove temporary files
echo $(date +"[%b %d %H:%M:%S] Removing temporary files...")

rm *txt.tmp

echo $(date +"[%b %d %H:%M:%S] Done")