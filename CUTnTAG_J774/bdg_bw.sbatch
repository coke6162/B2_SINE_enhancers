#!/bin/bash

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=48:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bdg_bw
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

chromSizesFile="/Shares/CL_Shared/db/genomes/mm10/fa/MusMus10_primary_chrom.sizes"
bedGraphFolder="/Users/%u/bdg"
outputFolder="/Users/%u/bigwig"

# first pass: sort the bedgraphs

for i in ${bedGraphFolder}/*treat*.bdg
do 
basename=`basename $i`
sort -k1,1 -k2,2n $i > ${outputFolder}/${basename/.bdg/.sorted.bdg} &
done
wait


# second pass: convert sorted bedgraphs to bigwigs

for i in ${outputFolder}/*.sorted.bdg
do
basename=`basename $i`
/Shares/CL_Shared/programs/ucsc/bedGraphToBigWig $i $chromSizesFile ${outputFolder}/${basename/.sorted.bdg/.bw} & 
done
wait

# OK to remove bedgraphs 
rm ${outputFolder}/*sorted.bdg
