#!/bin/bash

# Set seed for shuffling (adopted from https://stackoverflow.com/a/41962458/7820599)
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

seed=1;

# Option parsing adopted from https://stackoverflow.com/a/14203146
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

## Piccolo et al. 2017 IFNG 4h vs UT

# Get random 750 nonresponsive genes (padj > 0.90, baseMean > 100, abs(log2FC) < 0.10)
shuf -n 750 --random-source=<(get_seeded_random $seed) $REST picc_BMDM_IFNG_4h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.bed \
| bedtools sort -i \
> picc_BMDM_IFNG_4h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10_random750.bed

## Piccolo et al. 2017 IFNG 2h vs UT

# Get random 750 nonresponsive genes (padj > 0.90, baseMean > 100, abs(log2FC) < 0.10)
shuf -n 750 --random-source=<(get_seeded_random $seed) $REST picc_BMDM_IFNG_42_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.bed \
| bedtools sort -i \
> picc_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10_random750.bed

## Platanitis et al. 2019 IFNG 2h vs UT

# Get random 750 nonresponsive genes (padj > 0.90, baseMean > 100, abs(log2FC) < 0.10)
shuf -n 750 --random-source=<(get_seeded_random $seed) $REST plat_BMDM_IFNG_42_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.bed \
| bedtools sort -i \
> plat_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10_random750.bed