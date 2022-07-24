#!/bin/bash

# Write temporary file containing sample name as column 1 followed by GIGGLE output
for i in *.giggle.txt
do 
    cat ${i} | grep -v "#" | awk '{print "'$i'" "\t" $0}' | sed 's/_peaks.narrowPeak.giggle.txt//g' | sed 's/.giggle.txt//g' > ${i}.tmp
done

# Concatenate files
cat *.tmp \
| sort -nrk9 \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $9}' \
| sed '1isample\trepeat\tfilesize\toverlaps\toddsratio\tpvalue\tgigglescore' \
> all.giggle.txt

# Get names of families that meet minimum thresholds
awk '(NR>1){if ($3 > 100 && $4 > 30 && $5 > 3 && $7 > 100 && $1 !~ "UT") print $0}' all.giggle.txt \
| awk '{print $2}' \
> repeat_families_filtered.txt

# Identify subset families that are minimally enriched in H3K27ac (giggleScore > 1)
while read line
do
    awk -v var="$line" '{if ($2 == var && $1 ~ "picc_BMDM_WT_IFNG_2h_H3K27ac" && $7 > 1) print $2}' all.giggle.txt
done < repeat_families_filtered.txt | sort > repeat_families_filtered_H3K27ac.txt

# Identify subset families that are minimally enriched in STAT1 (giggleScore > 1)
while read line
do
    awk -v var="$line" '{if ($2 == var && $1 ~ "picc_BMDM_WT_IFNG_2h_STAT1" && $7 > 1) print $2}' all.giggle.txt
done < repeat_families_filtered.txt | sort > repeat_families_filtered_STAT1.txt

# Identify families enriched in both H3K27ac and STAT1
comm -12 repeat_families_filtered_H3K27ac.txt repeat_families_filtered_STAT1.txt
#B2_Mm2
#ID_B1
#L1MB2
#LTR104_Mam
#MIR1_Amn
#MIR3
#MT2B2
#PB1D7
#RLTR30B_MM
#RLTR30E_MM
#RLTR30D_MM

# Extract Giggle results for these families
cat all.giggle.txt \
| awk '{if ($2=="RLTR30B_MM" || $2=="B2_Mm2" || $2=="FAM" || $2=="ID_B1" || $2=="Kanga2_a" || $2=="L1MB2" || $2=="Looper" || $2=="LTR104_Mam" || $2=="LTR78" || $2=="MamRep605b" || $2=="MER102a" || $2=="MER102b" || $2=="MER45A" || $2=="MER96" || $2=="MER97a" || $2=="MIR1_Amn" || $2=="MIR3" || $2=="MT2B2" || $2=="PB1D7" || $2=="RLTR30D_MM" || $2=="RLTR30E_MM") print $0}' \
| sed 1i"sample\trepeat_family\trepeat_family_size\toverlaps\todds_ratio\tpvalue\tgiggle_score" \
| sed 's/_BMDM//g' \
> all_filtered.giggle.txt

# Remove temporary files
rm *.giggle.txt.tmp