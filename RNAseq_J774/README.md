The RNAseq workflow is as follows:
  1. bbduk_PE.sbatch
  2. fastqcreport.sbatch
  3. multiqc.sbatch
  4. hisat2_PE.sbatch
  5. merge_bams.sbatch
  6. bamTobw.sbatch
  7. featureCounts.sbatch
  8. deseq2.R

Distance plots graphing log2FC versus genome coordinates were used to visualize differential gene expression from DESeq2 output tables using the script distance_plots.R
